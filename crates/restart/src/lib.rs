//! Restart files: a run's state as named arrays of `f64`, saved and loaded
//! bit for bit.
//!
//! A restart holds state only, never settings: a restarted run builds its
//! world and model from the command line as usual, then each component
//! loads its arrays over the fresh state. The header carries a fingerprint
//! of the grid the arrays live on, so a file from another grid is refused,
//! and the number of steps taken when it was saved.
//!
//! # Format
//!
//! All little-endian:
//!
//! ```text
//! magic      8 bytes  b"FMSRSTRT"
//! version    u32      1
//! arrays     u32      number of arrays
//! grid       u64      grid fingerprint
//! step       u64      steps taken when saved
//! then per array, in insertion order:
//!   name length  u16
//!   name         UTF-8 bytes (non-empty, unique in the file)
//!   count        u64   number of values
//!   values       count × f64 (the bits, to_le_bytes)
//! ```
//!
//! Nothing may follow the last array. Values round-trip exactly, NaN
//! payloads and signed zeros included.

use std::collections::HashSet;
use std::fmt;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

/// The first eight bytes of every restart file.
pub const MAGIC: [u8; 8] = *b"FMSRSTRT";

/// The format version this build writes and reads.
pub const VERSION: u32 = 1;

/// A run's state: a header (grid fingerprint and step) and uniquely named
/// arrays of `f64`, kept in insertion order.
#[derive(Clone, Debug)]
pub struct Restart {
    grid: u64,
    step: u64,
    arrays: Vec<(String, Vec<f64>)>,
}

impl Restart {
    /// An empty restart for the grid with fingerprint `grid`, saved after
    /// `step` steps.
    pub fn new(grid: u64, step: u64) -> Self {
        Self {
            grid,
            step,
            arrays: Vec::new(),
        }
    }

    /// Fingerprint of the grid the arrays live on.
    pub fn grid(&self) -> u64 {
        self.grid
    }

    /// Steps taken when saved.
    pub fn step(&self) -> u64 {
        self.step
    }

    /// Add an array.
    ///
    /// # Panics
    /// If the name is empty, longer than `u16::MAX` bytes, or already used:
    /// these are programming errors.
    pub fn put(&mut self, name: &str, values: &[f64]) {
        assert!(!name.is_empty(), "a restart array needs a name");
        assert!(
            name.len() <= usize::from(u16::MAX),
            "restart array name longer than {} bytes",
            u16::MAX
        );
        assert!(
            !self.contains(name),
            "restart array `{name}` is already present"
        );
        self.arrays.push((name.to_owned(), values.to_vec()));
    }

    /// Whether there is an array called `name`.
    pub fn contains(&self, name: &str) -> bool {
        self.arrays.iter().any(|(n, _)| n == name)
    }

    /// The array called `name`, which must have `len` values.
    pub fn get(&self, name: &str, len: usize) -> Result<&[f64], RestartError> {
        let (_, values) =
            self.arrays
                .iter()
                .find(|(n, _)| n == name)
                .ok_or_else(|| RestartError::Missing {
                    name: name.to_owned(),
                })?;
        if values.len() == len {
            Ok(values)
        } else {
            Err(RestartError::WrongLength {
                name: name.to_owned(),
                expected: len,
                found: values.len(),
            })
        }
    }

    /// The array called `name`, which must have `N` values, as an array: for
    /// a short, fixed list of scalars.
    pub fn get_array<const N: usize>(&self, name: &str) -> Result<[f64; N], RestartError> {
        let values = self.get(name, N)?;
        Ok(std::array::from_fn(|i| values[i]))
    }

    /// The array called `name` as `rows` rows of `width` values each, which
    /// is how a list of fixed-width records is saved (row-major).
    ///
    /// # Panics
    /// If `width` is 0 or `rows × width` overflows: programming errors.
    pub fn get_rows(
        &self,
        name: &str,
        rows: usize,
        width: usize,
    ) -> Result<std::slice::ChunksExact<'_, f64>, RestartError> {
        assert!(width > 0, "a restart row needs a width");
        let len = rows
            .checked_mul(width)
            .expect("restart rows overflow a length");
        Ok(self.get(name, len)?.chunks_exact(width))
    }

    /// The arrays' names, in insertion order.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.arrays.iter().map(|(n, _)| n.as_str())
    }

    /// Every array with its name, in insertion order: for inspecting or
    /// copying a restart without knowing its arrays' lengths.
    pub fn arrays(&self) -> impl Iterator<Item = (&str, &[f64])> {
        self.arrays.iter().map(|(n, v)| (n.as_str(), v.as_slice()))
    }

    /// Write in the restart format (see the crate docs).
    pub fn write_to(&self, w: &mut impl Write) -> io::Result<()> {
        let count = u32::try_from(self.arrays.len())
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "too many restart arrays"))?;
        w.write_all(&MAGIC)?;
        w.write_all(&VERSION.to_le_bytes())?;
        w.write_all(&count.to_le_bytes())?;
        w.write_all(&self.grid.to_le_bytes())?;
        w.write_all(&self.step.to_le_bytes())?;
        let mut buffer = Vec::with_capacity(8 * CHUNK);
        for (name, values) in &self.arrays {
            // `put` keeps names within u16.
            w.write_all(&(name.len() as u16).to_le_bytes())?;
            w.write_all(name.as_bytes())?;
            w.write_all(&(values.len() as u64).to_le_bytes())?;
            for chunk in values.chunks(CHUNK) {
                buffer.clear();
                for v in chunk {
                    buffer.extend_from_slice(&v.to_bits().to_le_bytes());
                }
                w.write_all(&buffer)?;
            }
        }
        Ok(())
    }

    /// Read a restart written by [`write_to`](Self::write_to). The input
    /// must end where the restart does.
    pub fn read_from(r: &mut impl Read) -> Result<Self, RestartError> {
        let mut magic = Vec::with_capacity(MAGIC.len());
        r.by_ref()
            .take(MAGIC.len() as u64)
            .read_to_end(&mut magic)?;
        if magic.as_slice() != MAGIC {
            return Err(if MAGIC.starts_with(&magic) && !magic.is_empty() {
                truncated()
            } else {
                RestartError::NotRestart
            });
        }
        let version = u32::from_le_bytes(read_array(r)?);
        if version != VERSION {
            return Err(RestartError::UnsupportedVersion(version));
        }
        let count = u32::from_le_bytes(read_array(r)?);
        let grid = u64::from_le_bytes(read_array(r)?);
        let step = u64::from_le_bytes(read_array(r)?);
        let mut restart = Self::new(grid, step);
        let mut seen = HashSet::new();
        for _ in 0..count {
            let length = u16::from_le_bytes(read_array(r)?);
            if length == 0 {
                return Err(RestartError::Corrupt("an array has an empty name".into()));
            }
            let mut name = vec![0; usize::from(length)];
            read_exact(r, &mut name)?;
            let name = String::from_utf8(name)
                .map_err(|_| RestartError::Corrupt("an array name is not valid UTF-8".into()))?;
            if !seen.insert(name.clone()) {
                return Err(RestartError::Corrupt(format!(
                    "the array `{name}` appears twice"
                )));
            }
            let values = u64::from_le_bytes(read_array(r)?);
            let bytes = values.checked_mul(8).ok_or_else(|| {
                RestartError::Corrupt(format!("the array `{name}` has an impossible length"))
            })?;
            // Grows with what is actually read, never with the stored count.
            let mut raw = Vec::new();
            r.by_ref().take(bytes).read_to_end(&mut raw)?;
            if raw.len() as u64 != bytes {
                return Err(truncated());
            }
            let values = raw
                .as_chunks::<8>()
                .0
                .iter()
                .map(|&b| f64::from_bits(u64::from_le_bytes(b)))
                .collect();
            restart.arrays.push((name, values));
        }
        let mut extra = [0u8; 1];
        loop {
            match r.read(&mut extra) {
                Ok(0) => break,
                Ok(_) => {
                    return Err(RestartError::Corrupt(
                        "unexpected bytes after the last array".into(),
                    ));
                }
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e.into()),
            }
        }
        Ok(restart)
    }

    /// The restart in its file format.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        self.write_to(&mut bytes)
            .expect("writing to memory cannot fail");
        bytes
    }

    /// Save to `path`. Writes a temporary file beside it, then renames it
    /// over `path`, so a failed save never replaces a good file.
    pub fn save(&self, path: &Path) -> io::Result<()> {
        let name = path.file_name().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("not a file path: {}", path.display()),
            )
        })?;
        let mut temporary_name = name.to_owned();
        temporary_name.push(format!(".tmp-{}", std::process::id()));
        let temporary = path.with_file_name(temporary_name);
        let written = (|| {
            let mut w = BufWriter::new(File::create(&temporary)?);
            self.write_to(&mut w)?;
            w.into_inner()
                .map_err(io::IntoInnerError::into_error)?
                .sync_all()?;
            std::fs::rename(&temporary, path)
        })();
        if written.is_err() {
            // Best effort: the error that matters is the one returned.
            let _ = std::fs::remove_file(&temporary);
        }
        written
    }

    /// Load a restart saved by [`save`](Self::save).
    pub fn load(path: &Path) -> Result<Self, RestartError> {
        let mut r = io::BufReader::new(File::open(path)?);
        Self::read_from(&mut r)
    }
}

/// A counter saved as `f64` (in the array `name`), which must be a whole
/// number that fits `T`; otherwise [`RestartError::Mismatch`].
pub fn count<T: TryFrom<u64>>(name: &str, value: f64) -> Result<T, RestartError> {
    let not_a_count = || RestartError::Mismatch(format!("{name} holds {value}, not a count"));
    // 2⁶⁴ as f64: every whole number below it fits a u64.
    if (0.0..18_446_744_073_709_551_616.0).contains(&value) && value.fract() == 0.0 {
        T::try_from(value as u64).map_err(|_| not_a_count())
    } else {
        Err(not_a_count())
    }
}

/// A flag saved as `f64` (in the array `name`): 0 is false and 1 true;
/// anything else is [`RestartError::Mismatch`].
pub fn flag(name: &str, value: f64) -> Result<bool, RestartError> {
    if value == 0.0 {
        Ok(false)
    } else if value == 1.0 {
        Ok(true)
    } else {
        Err(RestartError::Mismatch(format!(
            "{name} holds {value}, not a flag (0 or 1)"
        )))
    }
}

/// Values per write.
const CHUNK: usize = 1024;

/// Read exactly `buf.len()` bytes; running out is corruption.
fn read_exact(r: &mut impl Read, buf: &mut [u8]) -> Result<(), RestartError> {
    r.read_exact(buf).map_err(|e| match e.kind() {
        io::ErrorKind::UnexpectedEof => truncated(),
        _ => e.into(),
    })
}

fn read_array<const N: usize>(r: &mut impl Read) -> Result<[u8; N], RestartError> {
    let mut buf = [0; N];
    read_exact(r, &mut buf)?;
    Ok(buf)
}

fn truncated() -> RestartError {
    RestartError::Corrupt("the file ends early (truncated)".into())
}

/// Why a restart could not be read or loaded.
#[derive(Debug)]
pub enum RestartError {
    /// Reading the file failed.
    Io(io::Error),
    /// The file does not start with [`MAGIC`].
    NotRestart,
    /// The file is in a format version this build doesn't read.
    UnsupportedVersion(u32),
    /// The file is damaged: truncated, an invalid or repeated name, or bytes
    /// after the last array.
    Corrupt(String),
    /// An array a component needs is not in the file.
    Missing {
        /// The array's name.
        name: String,
    },
    /// An array has the wrong number of values.
    WrongLength {
        /// The array's name.
        name: String,
        /// Values the run needs.
        expected: usize,
        /// Values in the file.
        found: usize,
    },
    /// The file was saved on another grid.
    DifferentGrid {
        /// The file's grid fingerprint.
        file: u64,
        /// This run's grid fingerprint.
        run: u64,
    },
    /// A component refused the file for its own reason (such as different
    /// vertical levels).
    Mismatch(String),
}

impl fmt::Display for RestartError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "restart file: {e}"),
            Self::NotRestart => write!(f, "not a restart file (it does not start with FMSRSTRT)"),
            Self::UnsupportedVersion(v) => write!(
                f,
                "restart file format version {v} is not supported (this build reads version {VERSION})"
            ),
            Self::Corrupt(why) => write!(f, "corrupt restart file: {why}"),
            Self::Missing { name } => write!(f, "the restart file has no array `{name}`"),
            Self::WrongLength {
                name,
                expected,
                found,
            } => write!(
                f,
                "the restart file's array `{name}` has {found} values; this run needs {expected}"
            ),
            Self::DifferentGrid { file, run } => write!(
                f,
                "the restart file was saved on a different grid (fingerprint {file:016x}; this run's is {run:016x})"
            ),
            Self::Mismatch(why) => write!(f, "the restart file does not fit this run: {why}"),
        }
    }
}

impl std::error::Error for RestartError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for RestartError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Restart {
        let mut r = Restart::new(0x0123_4567_89ab_cdef, 42);
        r.put(
            "b.first",
            &[1.0, -0.0, f64::NAN, f64::from_bits(0x7ff8_dead_beef_0001)],
        );
        r.put("a.second", &[]);
        r.put(
            "z.third",
            &[
                f64::MIN_POSITIVE / 3.0,
                f64::INFINITY,
                f64::NEG_INFINITY,
                -f64::from_bits(1),
            ],
        );
        r
    }

    fn bits(r: &Restart, name: &str, len: usize) -> Vec<u64> {
        r.get(name, len)
            .unwrap()
            .iter()
            .map(|v| v.to_bits())
            .collect()
    }

    #[test]
    fn round_trip_is_bit_exact_and_ordered() {
        let r = sample();
        let bytes = r.to_bytes();
        let back = Restart::read_from(&mut bytes.as_slice()).unwrap();
        assert_eq!(back.grid(), r.grid());
        assert_eq!(back.step(), 42);
        assert_eq!(
            back.names().collect::<Vec<_>>(),
            ["b.first", "a.second", "z.third"]
        );
        let lengths: Vec<(&str, usize)> = back.arrays().map(|(n, v)| (n, v.len())).collect();
        assert_eq!(lengths, [("b.first", 4), ("a.second", 0), ("z.third", 4)]);
        for (name, len) in [("b.first", 4), ("a.second", 0), ("z.third", 4)] {
            assert_eq!(bits(&back, name, len), bits(&r, name, len), "{name}");
        }
        assert_eq!(back.to_bytes(), bytes);
    }

    #[test]
    fn header_layout() {
        let bytes = Restart::new(7, 9).to_bytes();
        assert_eq!(bytes.len(), 32);
        assert_eq!(&bytes[..8], b"FMSRSTRT");
        assert_eq!(bytes[8..12], 1u32.to_le_bytes());
        assert_eq!(bytes[12..16], 0u32.to_le_bytes());
        assert_eq!(bytes[16..24], 7u64.to_le_bytes());
        assert_eq!(bytes[24..32], 9u64.to_le_bytes());
    }

    #[test]
    fn bad_magic_is_not_a_restart() {
        let mut bytes = sample().to_bytes();
        bytes[0] = b'X';
        assert!(matches!(
            Restart::read_from(&mut bytes.as_slice()),
            Err(RestartError::NotRestart)
        ));
        assert!(matches!(
            Restart::read_from(&mut &b""[..]),
            Err(RestartError::NotRestart)
        ));
    }

    #[test]
    fn other_versions_are_refused() {
        let mut bytes = sample().to_bytes();
        bytes[8..12].copy_from_slice(&2u32.to_le_bytes());
        assert!(matches!(
            Restart::read_from(&mut bytes.as_slice()),
            Err(RestartError::UnsupportedVersion(2))
        ));
    }

    #[test]
    fn every_truncation_fails_cleanly() {
        let bytes = sample().to_bytes();
        for end in 1..bytes.len() {
            match Restart::read_from(&mut &bytes[..end]) {
                Err(RestartError::Corrupt(_)) => {}
                other => panic!("prefix of {end} bytes: {other:?}"),
            }
        }
    }

    #[test]
    fn trailing_bytes_are_corrupt() {
        let mut bytes = sample().to_bytes();
        bytes.push(0);
        assert!(matches!(
            Restart::read_from(&mut bytes.as_slice()),
            Err(RestartError::Corrupt(_))
        ));
    }

    #[test]
    fn duplicate_names_in_a_file_are_corrupt() {
        let mut r = Restart::new(1, 2);
        r.put("x", &[1.0]);
        r.put("y", &[2.0]);
        let mut bytes = r.to_bytes();
        // Rename "y" to "x": the second array's name byte.
        let at = bytes.iter().rposition(|&b| b == b'y').unwrap();
        bytes[at] = b'x';
        let err = Restart::read_from(&mut bytes.as_slice()).unwrap_err();
        assert!(matches!(err, RestartError::Corrupt(_)), "{err:?}");
        assert!(err.to_string().contains("`x` appears twice"), "{err}");
    }

    #[test]
    fn invalid_names_are_corrupt() {
        let mut r = Restart::new(1, 2);
        r.put("x", &[1.0]);
        let mut bytes = r.to_bytes();
        let at = bytes.iter().position(|&b| b == b'x').unwrap();
        bytes[at] = 0xff;
        assert!(matches!(
            Restart::read_from(&mut bytes.as_slice()),
            Err(RestartError::Corrupt(_))
        ));
    }

    #[test]
    fn huge_counts_do_not_allocate() {
        let mut r = Restart::new(1, 2);
        r.put("x", &[1.0]);
        let mut bytes = r.to_bytes();
        // The value count, just before the last 8 bytes.
        let n = bytes.len();
        bytes[n - 16..n - 8].copy_from_slice(&(u64::MAX / 16).to_le_bytes());
        assert!(matches!(
            Restart::read_from(&mut bytes.as_slice()),
            Err(RestartError::Corrupt(_))
        ));
        bytes[n - 16..n - 8].copy_from_slice(&u64::MAX.to_le_bytes());
        assert!(matches!(
            Restart::read_from(&mut bytes.as_slice()),
            Err(RestartError::Corrupt(_))
        ));
        // Many arrays claimed, none present.
        let mut bytes = Restart::new(1, 2).to_bytes();
        bytes[12..16].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(matches!(
            Restart::read_from(&mut bytes.as_slice()),
            Err(RestartError::Corrupt(_))
        ));
    }

    #[test]
    fn get_checks_presence_and_length() {
        let r = sample();
        assert!(r.contains("b.first"));
        assert!(!r.contains("nope"));
        assert!(matches!(
            r.get("nope", 4),
            Err(RestartError::Missing { name }) if name == "nope"
        ));
        assert!(matches!(
            r.get("b.first", 5),
            Err(RestartError::WrongLength {
                expected: 5,
                found: 4,
                ..
            })
        ));
    }

    #[test]
    fn fixed_arrays_and_rows_check_their_length() {
        let r = sample();
        let a: [f64; 4] = r.get_array("z.third").unwrap();
        assert_eq!(a.map(f64::to_bits), bits(&r, "z.third", 4)[..]);
        assert!(matches!(
            r.get_array::<3>("z.third"),
            Err(RestartError::WrongLength {
                expected: 3,
                found: 4,
                ..
            })
        ));
        assert!(matches!(
            r.get_array::<1>("nope"),
            Err(RestartError::Missing { .. })
        ));
        let rows: Vec<&[f64]> = r.get_rows("z.third", 2, 2).unwrap().collect();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[1][0], f64::NEG_INFINITY);
        assert_eq!(r.get_rows("a.second", 0, 7).unwrap().count(), 0);
        assert!(matches!(
            r.get_rows("z.third", 3, 2),
            Err(RestartError::WrongLength {
                expected: 6,
                found: 4,
                ..
            })
        ));
    }

    #[test]
    fn counts_are_whole_numbers_in_range() {
        assert_eq!(count::<u64>("c", 0.0).unwrap(), 0);
        assert_eq!(count::<u32>("c", 4_294_967_295.0).unwrap(), u32::MAX);
        assert_eq!(count::<usize>("c", 12.0).unwrap(), 12);
        assert_eq!(count::<u64>("c", 9_007_199_254_740_992.0).unwrap(), 1 << 53);
        for bad in [
            -1.0,
            0.5,
            f64::NAN,
            f64::INFINITY,
            18_446_744_073_709_551_616.0,
        ] {
            let err = count::<u64>("c", bad).unwrap_err();
            assert!(matches!(err, RestartError::Mismatch(_)), "{bad}: {err:?}");
        }
        let err = count::<u32>("x.counts", 4_294_967_296.0).unwrap_err();
        assert!(err.to_string().contains("x.counts"), "{err}");
        // Negative zero is zero.
        assert_eq!(count::<u8>("c", -0.0).unwrap(), 0);
        assert!(count::<u8>("c", 256.0).is_err());
    }

    #[test]
    fn flags_are_zero_or_one() {
        assert!(!flag("f", 0.0).unwrap());
        assert!(!flag("f", -0.0).unwrap());
        assert!(flag("f", 1.0).unwrap());
        for bad in [2.0, 0.5, -1.0, f64::NAN] {
            assert!(matches!(flag("f", bad), Err(RestartError::Mismatch(_))));
        }
    }

    #[test]
    #[should_panic(expected = "already present")]
    fn put_refuses_duplicates() {
        let mut r = Restart::new(0, 0);
        r.put("x", &[1.0]);
        r.put("x", &[2.0]);
    }

    #[test]
    #[should_panic(expected = "needs a name")]
    fn put_refuses_empty_names() {
        Restart::new(0, 0).put("", &[1.0]);
    }

    #[test]
    fn messages_read_well() {
        let e = RestartError::DifferentGrid { file: 1, run: 2 };
        assert!(e.to_string().contains("different grid"), "{e}");
        let e = RestartError::from(io::Error::new(io::ErrorKind::NotFound, "gone"));
        assert!(std::error::Error::source(&e).is_some());
    }

    #[test]
    fn save_and_load_replace_atomically() {
        let dir = std::env::temp_dir();
        let name = format!(
            "fms-restart-test-{}-{}.fmsrst",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        );
        let path = dir.join(&name);
        let leftovers = || {
            std::fs::read_dir(&dir)
                .unwrap()
                .filter_map(|e| e.ok())
                .filter(|e| {
                    let n = e.file_name();
                    let n = n.to_string_lossy();
                    n.starts_with(&name) && *n != name
                })
                .count()
        };

        let first = sample();
        first.save(&path).unwrap();
        assert_eq!(Restart::load(&path).unwrap().to_bytes(), first.to_bytes());

        let mut second = Restart::new(5, 6);
        second.put("only", &[3.0]);
        second.save(&path).unwrap();
        let back = Restart::load(&path).unwrap();
        assert_eq!(back.to_bytes(), second.to_bytes());
        assert_eq!(leftovers(), 0);

        std::fs::remove_file(&path).unwrap();
        assert!(matches!(Restart::load(&path), Err(RestartError::Io(_))));
    }
}
