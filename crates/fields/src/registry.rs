//! The field registry: every per-element layer of world state, by name.

use std::collections::HashMap;
use std::fmt;
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};
use std::sync::{RwLock, RwLockReadGuard, RwLockWriteGuard, TryLockError};

use crate::column::{Column, DType, FieldType};
use crate::meta::{Buffering, FieldMeta, FieldSpec, GridSizes};

/// Untyped handle to a registered field.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct FieldKey(u32);

impl FieldKey {
    /// Registration index (stable for the registry's lifetime).
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Typed handle to a registered field of element type `T`.
///
/// Obtained from [`FieldRegistry::register`] or [`FieldRegistry::id`], which
/// check the type once, so accesses through it need no further checks.
pub struct FieldId<T> {
    key: FieldKey,
    _type: PhantomData<fn() -> T>,
}

impl<T> FieldId<T> {
    /// The untyped key.
    pub fn key(self) -> FieldKey {
        self.key
    }
}

impl<T> Clone for FieldId<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for FieldId<T> {}
impl<T> PartialEq for FieldId<T> {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
    }
}
impl<T> Eq for FieldId<T> {}
impl<T> fmt::Debug for FieldId<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FieldId({})", self.key.0)
    }
}

/// Errors from registering or looking up fields.
#[derive(Clone, Debug, PartialEq)]
pub enum FieldError {
    /// A field with this name is already registered.
    DuplicateName(&'static str),
    /// The name does not follow the `namespace.field` convention
    /// (see [`validate_name`]).
    InvalidName(&'static str),
    /// No field has this name.
    UnknownName(String),
    /// The field exists but has a different element type.
    TypeMismatch {
        /// Field name.
        name: String,
        /// Type it was registered with.
        registered: DType,
        /// Type that was requested.
        requested: DType,
    },
    /// The initial values do not match the field's location.
    WrongLength {
        /// Field name.
        name: &'static str,
        /// Length required by the location.
        expected: usize,
        /// Length supplied.
        actual: usize,
    },
}

impl fmt::Display for FieldError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FieldError::DuplicateName(n) => write!(f, "field `{n}` is already registered"),
            FieldError::InvalidName(n) => write!(
                f,
                "field name `{n}` must be `namespace.field`, both parts snake_case \
                 (e.g. `terrain.elevation`)"
            ),
            FieldError::UnknownName(n) => write!(f, "no field named `{n}`"),
            FieldError::TypeMismatch {
                name,
                registered,
                requested,
            } => write!(
                f,
                "field `{name}` has type {registered:?}, not {requested:?}"
            ),
            FieldError::WrongLength {
                name,
                expected,
                actual,
            } => write!(f, "field `{name}` needs {expected} values, got {actual}"),
        }
    }
}

impl std::error::Error for FieldError {}

struct Entry {
    meta: FieldMeta,
    front: RwLock<Column>,
    back: Option<RwLock<Column>>,
}

/// Every field of the world, with metadata and per-field locks.
///
/// Fields are registered up front (`&mut self`); after that, any number of threads
/// can access *different* fields concurrently through `&self`. The scheduler
/// guarantees that concurrently running systems never touch the same field
/// mutably. The per-field locks enforce this: a conflicting access **panics**
/// instead of blocking, because it means a system under-declared its reads or writes.
pub struct FieldRegistry {
    sizes: GridSizes,
    entries: Vec<Entry>,
    by_name: HashMap<&'static str, FieldKey>,
}

impl FieldRegistry {
    /// An empty registry for grids of the given sizes.
    pub fn new(sizes: GridSizes) -> Self {
        Self {
            sizes,
            entries: Vec::new(),
            by_name: HashMap::new(),
        }
    }

    /// Grid sizes fields are allocated against.
    pub fn sizes(&self) -> GridSizes {
        self.sizes
    }

    /// Register a field filled with `fill`.
    pub fn register<T: FieldType>(
        &mut self,
        spec: FieldSpec,
        fill: T,
    ) -> Result<FieldId<T>, FieldError> {
        let len = self.sizes.len(spec.location);
        self.register_with(spec, vec![fill; len])
    }

    /// Register a field with explicit initial values (copied to both buffers if
    /// double buffered).
    pub fn register_with<T: FieldType>(
        &mut self,
        spec: FieldSpec,
        values: Vec<T>,
    ) -> Result<FieldId<T>, FieldError> {
        if !validate_name(spec.name) {
            return Err(FieldError::InvalidName(spec.name));
        }
        if self.by_name.contains_key(spec.name) {
            return Err(FieldError::DuplicateName(spec.name));
        }
        let expected = self.sizes.len(spec.location);
        if values.len() != expected {
            return Err(FieldError::WrongLength {
                name: spec.name,
                expected,
                actual: values.len(),
            });
        }
        let back = match spec.buffering {
            Buffering::Single => None,
            Buffering::Double => Some(RwLock::new(T::into_column(values.clone()))),
        };
        let key = FieldKey(self.entries.len() as u32);
        self.by_name.insert(spec.name, key);
        self.entries.push(Entry {
            meta: FieldMeta {
                spec,
                dtype: T::DTYPE,
                len: expected,
            },
            front: RwLock::new(T::into_column(values)),
            back,
        });
        Ok(FieldId {
            key,
            _type: PhantomData,
        })
    }

    /// Typed handle for the field called `name`.
    pub fn id<T: FieldType>(&self, name: &str) -> Result<FieldId<T>, FieldError> {
        let key = self.key(name)?;
        let registered = self.entries[key.index()].meta.dtype;
        if registered != T::DTYPE {
            return Err(FieldError::TypeMismatch {
                name: name.to_owned(),
                registered,
                requested: T::DTYPE,
            });
        }
        Ok(FieldId {
            key,
            _type: PhantomData,
        })
    }

    /// Untyped key for the field called `name`.
    pub fn key(&self, name: &str) -> Result<FieldKey, FieldError> {
        self.by_name
            .get(name)
            .copied()
            .ok_or_else(|| FieldError::UnknownName(name.to_owned()))
    }

    /// Metadata of a field.
    pub fn meta(&self, key: FieldKey) -> &FieldMeta {
        &self.entries[key.index()].meta
    }

    /// All registered fields, in registration order.
    pub fn keys(&self) -> impl Iterator<Item = FieldKey> + '_ {
        (0..self.entries.len() as u32).map(FieldKey)
    }

    /// Number of registered fields.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether no fields are registered.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Read the (front buffer of the) field.
    ///
    /// # Panics
    /// If the field is currently borrowed mutably.
    pub fn read<T: FieldType>(&self, id: FieldId<T>) -> FieldRef<'_, T> {
        let entry = &self.entries[id.key.index()];
        FieldRef {
            guard: lock_read(&entry.front, &entry.meta),
            _type: PhantomData,
        }
    }

    /// Write the (front buffer of the) field in place.
    ///
    /// # Panics
    /// If the field is currently borrowed.
    pub fn write<T: FieldType>(&self, id: FieldId<T>) -> FieldMut<'_, T> {
        let entry = &self.entries[id.key.index()];
        FieldMut {
            guard: lock_write(&entry.front, &entry.meta),
            _type: PhantomData,
        }
    }

    /// Read the back buffer of a double-buffered field.
    ///
    /// # Panics
    /// If the field is single buffered, or its back buffer is borrowed mutably.
    pub fn read_back<T: FieldType>(&self, id: FieldId<T>) -> FieldRef<'_, T> {
        let entry = &self.entries[id.key.index()];
        FieldRef {
            guard: lock_read(back_buffer(entry), &entry.meta),
            _type: PhantomData,
        }
    }

    /// Write the back buffer of a double-buffered field.
    ///
    /// # Panics
    /// If the field is single buffered, or its back buffer is borrowed.
    pub fn write_back<T: FieldType>(&self, id: FieldId<T>) -> FieldMut<'_, T> {
        let entry = &self.entries[id.key.index()];
        FieldMut {
            guard: lock_write(back_buffer(entry), &entry.meta),
            _type: PhantomData,
        }
    }

    /// Exchange the front and back buffers of a double-buffered field. O(1).
    ///
    /// # Panics
    /// If the field is single buffered, or either buffer is borrowed.
    pub fn swap(&self, key: FieldKey) {
        let entry = &self.entries[key.index()];
        let mut back = lock_write(back_buffer(entry), &entry.meta);
        let mut front = lock_write(&entry.front, &entry.meta);
        std::mem::swap(&mut *front, &mut *back);
    }

    /// Untyped read access, e.g. for map modes that display any field.
    ///
    /// # Panics
    /// If the field is currently borrowed mutably.
    pub fn read_column(&self, key: FieldKey) -> ColumnRef<'_> {
        let entry = &self.entries[key.index()];
        ColumnRef {
            guard: lock_read(&entry.front, &entry.meta),
        }
    }

    /// First element that is non-finite or outside the field's valid range, as
    /// `(index, value)`. Used to detect numerical blow-ups early.
    pub fn first_invalid(&self, key: FieldKey) -> Option<(usize, f64)> {
        let range = self.meta(key).spec.valid_range;
        let column = self.read_column(key);
        (0..column.len())
            .map(|i| (i, column.get_f64(i)))
            .find(|&(_, v)| !v.is_finite() || range.is_some_and(|(lo, hi)| v < lo || v > hi))
    }
}

impl fmt::Debug for FieldRegistry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FieldRegistry")
            .field("sizes", &self.sizes)
            .field(
                "fields",
                &self
                    .entries
                    .iter()
                    .map(|e| e.meta.spec.name)
                    .collect::<Vec<_>>(),
            )
            .finish()
    }
}

/// Whether `name` follows the field naming convention: `namespace.field`, where
/// both parts are non-empty `snake_case` (lowercase ASCII letters, digits and
/// underscores, starting with a letter). The namespace names the owning system or
/// crate (`terrain`, `atmos`, `ocean`, `hydro`, `climate`, `civ`, …), so
/// independently written crates cannot collide.
pub fn validate_name(name: &str) -> bool {
    let part_ok = |p: &str| {
        let mut chars = p.chars();
        chars.next().is_some_and(|c| c.is_ascii_lowercase())
            && chars.all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_')
    };
    match name.split_once('.') {
        Some((namespace, field)) => part_ok(namespace) && part_ok(field),
        None => false,
    }
}

fn back_buffer(entry: &Entry) -> &RwLock<Column> {
    entry
        .back
        .as_ref()
        .unwrap_or_else(|| panic!("field `{}` is single buffered", entry.meta.spec.name))
}

fn lock_read<'a>(lock: &'a RwLock<Column>, meta: &FieldMeta) -> RwLockReadGuard<'a, Column> {
    match lock.try_read() {
        Ok(guard) => guard,
        Err(TryLockError::WouldBlock) => panic!(
            "access conflict: field `{}` is being written; a system's read/write set is \
             under-declared",
            meta.spec.name
        ),
        Err(TryLockError::Poisoned(_)) => {
            panic!(
                "field `{}` was poisoned by a panic during a write",
                meta.spec.name
            )
        }
    }
}

fn lock_write<'a>(lock: &'a RwLock<Column>, meta: &FieldMeta) -> RwLockWriteGuard<'a, Column> {
    match lock.try_write() {
        Ok(guard) => guard,
        Err(TryLockError::WouldBlock) => panic!(
            "access conflict: field `{}` is already borrowed; a system's read/write set is \
             under-declared",
            meta.spec.name
        ),
        Err(TryLockError::Poisoned(_)) => {
            panic!(
                "field `{}` was poisoned by a panic during a write",
                meta.spec.name
            )
        }
    }
}

/// Shared access to a field's values; derefs to `&[T]`.
pub struct FieldRef<'a, T: FieldType> {
    guard: RwLockReadGuard<'a, Column>,
    _type: PhantomData<T>,
}

impl<T: FieldType> Deref for FieldRef<'_, T> {
    type Target = [T];
    fn deref(&self) -> &[T] {
        T::slice(&self.guard).expect("type checked when the FieldId was created")
    }
}

/// Exclusive access to a field's values; derefs to `&mut [T]`.
pub struct FieldMut<'a, T: FieldType> {
    guard: RwLockWriteGuard<'a, Column>,
    _type: PhantomData<T>,
}

impl<T: FieldType> Deref for FieldMut<'_, T> {
    type Target = [T];
    fn deref(&self) -> &[T] {
        T::slice(&self.guard).expect("type checked when the FieldId was created")
    }
}

impl<T: FieldType> DerefMut for FieldMut<'_, T> {
    fn deref_mut(&mut self) -> &mut [T] {
        T::slice_mut(&mut self.guard).expect("type checked when the FieldId was created")
    }
}

/// Untyped shared access to a field; derefs to [`Column`].
pub struct ColumnRef<'a> {
    guard: RwLockReadGuard<'a, Column>,
}

impl Deref for ColumnRef<'_> {
    type Target = Column;
    fn deref(&self) -> &Column {
        &self.guard
    }
}
