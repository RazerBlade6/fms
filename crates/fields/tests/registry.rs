//! Field registry behaviour: registration, typed lookup, borrowing, buffers.

use fms_fields::{
    Buffering, Column, DType, FieldError, FieldRegistry, FieldSpec, GridSizes, Location, MeshCounts,
};
use rayon::prelude::*;

fn registry() -> FieldRegistry {
    FieldRegistry::new(GridSizes {
        terrain: MeshCounts {
            cells: 1000,
            edges: 2994,
            vertices: 1996,
        },
        climate: MeshCounts {
            cells: 42,
            edges: 120,
            vertices: 80,
        },
    })
}

#[test]
fn fields_have_the_length_of_their_location() {
    let mut r = registry();
    let a = r
        .register(
            FieldSpec::new("elevation", "m", Location::TERRAIN_CELLS),
            0.0f32,
        )
        .unwrap();
    let b = r
        .register(
            FieldSpec::new("u", "m s-1", Location::CLIMATE_EDGES),
            0.0f64,
        )
        .unwrap();
    let c = r
        .register(
            FieldSpec::new("zeta", "s-1", Location::CLIMATE_VERTICES),
            0.0f64,
        )
        .unwrap();
    assert_eq!(r.read(a).len(), 1000);
    assert_eq!(r.read(b).len(), 120);
    assert_eq!(r.read(c).len(), 80);
    assert_eq!(r.meta(a.key()).len, 1000);
    assert_eq!(r.meta(a.key()).dtype, DType::F32);
}

#[test]
fn names_are_unique_and_lookups_are_type_checked() {
    let mut r = registry();
    let spec = FieldSpec::new("koppen", "", Location::TERRAIN_CELLS);
    let id = r.register(spec.clone(), 0u8).unwrap();
    assert_eq!(
        r.register(spec, 0u8),
        Err(FieldError::DuplicateName("koppen"))
    );
    assert_eq!(r.id::<u8>("koppen"), Ok(id));
    assert!(matches!(
        r.id::<f32>("koppen"),
        Err(FieldError::TypeMismatch { .. })
    ));
    assert!(matches!(
        r.id::<u8>("nope"),
        Err(FieldError::UnknownName(_))
    ));
    assert_eq!(r.len(), 1);
}

#[test]
fn register_with_checks_length() {
    let mut r = registry();
    let err = r
        .register_with(
            FieldSpec::new("x", "", Location::CLIMATE_CELLS),
            vec![1.0f64; 7],
        )
        .unwrap_err();
    assert_eq!(
        err,
        FieldError::WrongLength {
            name: "x",
            expected: 42,
            actual: 7
        }
    );
    let id = r
        .register_with(
            FieldSpec::new("x", "", Location::CLIMATE_CELLS),
            vec![2.5f64; 42],
        )
        .unwrap();
    assert!(r.read(id).iter().all(|&v| v == 2.5));
}

#[test]
fn different_fields_can_be_accessed_concurrently() {
    let mut r = registry();
    let ids: Vec<_> = (0..8)
        .map(|i| {
            let name: &'static str = Box::leak(format!("field_{i}").into_boxed_str());
            r.register(FieldSpec::new(name, "", Location::TERRAIN_CELLS), 0.0f64)
                .unwrap()
        })
        .collect();
    let source = r
        .register(
            FieldSpec::new("source", "", Location::TERRAIN_CELLS),
            1.0f64,
        )
        .unwrap();
    // Eight writers in parallel, all reading the same shared source field.
    ids.par_iter().enumerate().for_each(|(i, &id)| {
        let src = r.read(source);
        let mut out = r.write(id);
        for (o, s) in out.iter_mut().zip(src.iter()) {
            *o = s * i as f64;
        }
    });
    for (i, &id) in ids.iter().enumerate() {
        assert!(r.read(id).iter().all(|&v| v == i as f64));
    }
}

#[test]
#[should_panic(expected = "access conflict: field `t` is already borrowed")]
fn conflicting_access_panics_instead_of_blocking() {
    let mut r = registry();
    let t = r
        .register(FieldSpec::new("t", "K", Location::CLIMATE_CELLS), 0.0f64)
        .unwrap();
    let _reader = r.read(t);
    let _writer = r.write(t);
}

#[test]
fn double_buffers_swap() {
    let mut r = registry();
    let spec = FieldSpec {
        buffering: Buffering::Double,
        ..FieldSpec::new("h", "m", Location::CLIMATE_CELLS)
    };
    let h = r.register(spec, 1.0f64).unwrap();
    {
        let current = r.read(h);
        let mut next = r.write_back(h);
        for (n, c) in next.iter_mut().zip(current.iter()) {
            *n = c + 1.0;
        }
    }
    assert!(
        r.read(h).iter().all(|&v| v == 1.0),
        "front unchanged before swap"
    );
    r.swap(h.key());
    assert!(r.read(h).iter().all(|&v| v == 2.0));
    assert!(r.read_back(h).iter().all(|&v| v == 1.0));
}

#[test]
#[should_panic(expected = "single buffered")]
fn back_buffer_of_single_buffered_field_panics() {
    let mut r = registry();
    let t = r
        .register(FieldSpec::new("t", "K", Location::CLIMATE_CELLS), 0.0f64)
        .unwrap();
    let _ = r.write_back(t);
}

#[test]
fn untyped_access_and_validity_checks() {
    let mut r = registry();
    let spec = FieldSpec {
        valid_range: Some((150.0, 350.0)),
        ..FieldSpec::new("t", "K", Location::CLIMATE_CELLS)
    };
    let t = r.register(spec, 288.0f32).unwrap();
    let k = r
        .register(FieldSpec::new("class", "", Location::CLIMATE_CELLS), 3u8)
        .unwrap();

    assert_eq!(r.read_column(k.key()).get_f64(5), 3.0);
    assert!(matches!(&*r.read_column(t.key()), Column::F32(v) if v.len() == 42));
    assert_eq!(r.first_invalid(t.key()), None);

    r.write(t)[10] = 400.0;
    assert_eq!(r.first_invalid(t.key()), Some((10, 400.0)));
    r.write(t)[10] = 300.0;
    r.write(t)[20] = f32::NAN;
    assert!(matches!(r.first_invalid(t.key()), Some((20, v)) if v.is_nan()));
}

#[test]
fn keys_iterate_in_registration_order() {
    let mut r = registry();
    for name in ["a", "b", "c"] {
        r.register(FieldSpec::new(name, "", Location::CLIMATE_CELLS), 0i32)
            .unwrap();
    }
    let names: Vec<_> = r.keys().map(|k| r.meta(k).spec.name).collect();
    assert_eq!(names, ["a", "b", "c"]);
}
