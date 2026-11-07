//! Typed per-element storage.

/// Element type of a field.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum DType {
    /// 32-bit float.
    F32,
    /// 64-bit float.
    F64,
    /// Unsigned byte (e.g. a class index such as a Köppen code).
    U8,
    /// Unsigned 16-bit integer.
    U16,
    /// Unsigned 32-bit integer (e.g. a cell or plate index).
    U32,
    /// Signed 32-bit integer.
    I32,
}

/// The values of one field: one contiguous array (struct-of-arrays layout).
#[derive(Clone, Debug, PartialEq)]
pub enum Column {
    /// 32-bit float values.
    F32(Vec<f32>),
    /// 64-bit float values.
    F64(Vec<f64>),
    /// Unsigned byte values.
    U8(Vec<u8>),
    /// Unsigned 16-bit values.
    U16(Vec<u16>),
    /// Unsigned 32-bit values.
    U32(Vec<u32>),
    /// Signed 32-bit values.
    I32(Vec<i32>),
}

impl Column {
    /// Element type.
    pub fn dtype(&self) -> DType {
        match self {
            Column::F32(_) => DType::F32,
            Column::F64(_) => DType::F64,
            Column::U8(_) => DType::U8,
            Column::U16(_) => DType::U16,
            Column::U32(_) => DType::U32,
            Column::I32(_) => DType::I32,
        }
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        match self {
            Column::F32(v) => v.len(),
            Column::F64(v) => v.len(),
            Column::U8(v) => v.len(),
            Column::U16(v) => v.len(),
            Column::U32(v) => v.len(),
            Column::I32(v) => v.len(),
        }
    }

    /// Whether the column has no elements.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Element `i` converted to `f64`, for type-agnostic consumers such as map
    /// modes and diagnostics.
    pub fn get_f64(&self, i: usize) -> f64 {
        match self {
            Column::F32(v) => f64::from(v[i]),
            Column::F64(v) => v[i],
            Column::U8(v) => f64::from(v[i]),
            Column::U16(v) => f64::from(v[i]),
            Column::U32(v) => f64::from(v[i]),
            Column::I32(v) => f64::from(v[i]),
        }
    }

    /// All elements converted to `f64`.
    pub fn to_f64_vec(&self) -> Vec<f64> {
        (0..self.len()).map(|i| self.get_f64(i)).collect()
    }
}

/// Rust types that can be stored in a field.
pub trait FieldType: Copy + Send + Sync + 'static {
    /// The matching [`DType`].
    const DTYPE: DType;
    /// Wrap a vector into a [`Column`].
    fn into_column(values: Vec<Self>) -> Column;
    /// View a column as a slice of `Self`, if the types match.
    fn slice(column: &Column) -> Option<&[Self]>;
    /// Mutable view of a column, if the types match.
    fn slice_mut(column: &mut Column) -> Option<&mut [Self]>;
}

macro_rules! field_type {
    ($t:ty, $variant:ident) => {
        impl FieldType for $t {
            const DTYPE: DType = DType::$variant;
            fn into_column(values: Vec<Self>) -> Column {
                Column::$variant(values)
            }
            fn slice(column: &Column) -> Option<&[Self]> {
                match column {
                    Column::$variant(v) => Some(v),
                    _ => None,
                }
            }
            fn slice_mut(column: &mut Column) -> Option<&mut [Self]> {
                match column {
                    Column::$variant(v) => Some(v),
                    _ => None,
                }
            }
        }
    };
}

field_type!(f32, F32);
field_type!(f64, F64);
field_type!(u8, U8);
field_type!(u16, U16);
field_type!(u32, U32);
field_type!(i32, I32);
