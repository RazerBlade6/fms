//! `fms-climate climate --save-restart` and `--restart`, run as the real
//! binary on the smallest world (162 cells, a 20-day year): a run resumed
//! from a restart ends where the straight run does, byte for byte, across
//! an ocean jump too; runs with other settings start warm and say so; files
//! that don't fit are refused.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

/// A fresh directory under the system's temporary one, removed when
/// dropped.
struct TempDir(PathBuf);

impl TempDir {
    fn new(name: &str) -> Self {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "fms-climate-restart-{name}-{}-{nanos}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    fn join(&self, path: &str) -> PathBuf {
        self.0.join(path)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The smallest world and output, but for the climate grid.
const SMALL: [&str; 10] = [
    "--terrain-frequency",
    "8",
    "--year-days",
    "20",
    "--width",
    "32",
    "--frames",
    "2",
    "--movie-days",
    "2",
];

/// `fms-climate climate` with the small world on a climate grid of
/// `frequency` (162 cells at 4), `args` and its output in `out`: whether it
/// succeeded, and its standard error.
fn climate(frequency: &str, out: &Path, args: &[&str]) -> (bool, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_fms-climate"))
        .arg("climate")
        .args(SMALL)
        .args(["--climate-frequency", frequency])
        .args(args)
        .arg("--out-dir")
        .arg(out)
        .output()
        .expect("the binary runs");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    (output.status.success(), stderr)
}

/// As [`climate`], which must succeed; returns its standard error.
fn run(out: &Path, args: &[&str]) -> String {
    let (ok, stderr) = climate("4", out, args);
    assert!(ok, "fms-climate climate {args:?} failed:\n{stderr}");
    stderr
}

fn read(path: &Path) -> Vec<u8> {
    std::fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn text(path: &Path) -> String {
    String::from_utf8(read(path)).unwrap()
}

/// The two runs' outputs are the same byte for byte.
fn assert_same(a: &Path, b: &Path) {
    assert!(
        read(a) == read(b),
        "{} and {} differ",
        a.display(),
        b.display()
    );
}

fn arg(path: &Path) -> &str {
    path.to_str().expect("a UTF-8 temporary path")
}

#[test]
fn a_resumed_run_continues_exactly_and_others_start_warm() {
    let dir = TempDir::new("exact");
    let (straight, part, resumed) = (dir.join("straight"), dir.join("part"), dir.join("resumed"));
    let true_time = ["--ocean-jumps", "0"];
    let saved = |name: &str| dir.join(&format!("{name}.fmsrst"));
    let (straight_file, part_file, resumed_file) =
        (saved("straight"), saved("part"), saved("resumed"));
    std::thread::scope(|s| {
        s.spawn(|| {
            let stderr = run(
                &straight,
                &[
                    &true_time[..],
                    &["--years", "3", "--save-restart", arg(&straight_file)],
                ]
                .concat(),
            );
            assert!(!stderr.contains("restart:"), "{stderr}");
        });
        run(
            &part,
            &[
                &true_time[..],
                &["--years", "2", "--save-restart", arg(&part_file)],
            ]
            .concat(),
        );
    });
    let stderr = run(
        &resumed,
        &[
            &true_time[..],
            &[
                "--restart",
                arg(&part_file),
                "--years",
                "1",
                "--save-restart",
                arg(&resumed_file),
            ],
        ]
        .concat(),
    );
    assert!(
        stderr.contains(&format!(
            "restart: an exact continuation of {} from step 960 (day 40.0, in year 3)",
            part_file.display()
        )),
        "{stderr}"
    );
    assert_same(&straight_file, &resumed_file);
    for name in ["annual_cells.csv", "zonal_mean.csv", "overturning.csv"] {
        assert_same(
            &straight.join("climate").join(name),
            &resumed.join("climate").join(name),
        );
    }
    let summary = text(&resumed.join("climate/summary.txt"));
    assert!(
        summary.contains("restart: an exact continuation"),
        "{summary}"
    );

    // The monthly means go on from the third year's first month: the
    // straight run's last year, row for row.
    let monthly = |out: &Path| text(&out.join("climate/monthly_zonal.csv"));
    let (all, last_year) = (monthly(&straight), monthly(&resumed));
    let header = all.lines().next().unwrap();
    let year_3: Vec<&str> = all.lines().filter(|l| l.starts_with("3,")).collect();
    let resumed_rows: Vec<&str> = last_year.lines().collect();
    assert_eq!(resumed_rows[0], header);
    assert!(resumed_rows[1].starts_with("3,1,"), "{}", resumed_rows[1]);
    assert!(resumed_rows.last().unwrap().starts_with("3,12,"));
    assert!(resumed_rows[1..] == year_3[..], "the resumed months differ");

    // Warm starts from the same file: another CO₂, no gyres, another step.
    std::thread::scope(|s| {
        for (name, options, note) in [
            (
                "co2",
                &["--co2", "840"][..],
                "settings differ from the saved run's",
            ),
            (
                "gyres",
                &["--no-ocean-gyres"][..],
                "the wind-driven gyres are off in this run and were on in the saved one",
            ),
            (
                "step",
                &["--max-wind", "300"][..],
                "the time step changed from 3600 s to 3456 s",
            ),
        ] {
            let (dir, part_file) = (&dir, &part_file);
            s.spawn(move || {
                let args = [
                    &true_time[..],
                    &["--restart", arg(part_file), "--years", "1"],
                    options,
                ]
                .concat();
                let stderr = run(&dir.join(name), &args);
                assert!(
                    stderr.contains("not an exact continuation") && stderr.contains(note),
                    "{name}: {stderr}"
                );
                eprintln!(
                    "{name}: {}",
                    stderr
                        .lines()
                        .skip_while(|l| !l.contains("restart:"))
                        .take_while(|l| !l.starts_with("ran in"))
                        .collect::<Vec<_>>()
                        .join("\n")
                );
            });
        }
    });
}

#[test]
fn a_run_resumed_before_an_ocean_jump_continues_exactly() {
    // Land and coasts, and a jump at the end of year 2: saved half way
    // through year 2, which is being measured for it.
    let dir = TempDir::new("jump");
    let schedule = [
        "--continent",
        "--ocean-jumps",
        "1",
        "--jump-first-year",
        "2",
    ];
    let saved = |name: &str| dir.join(&format!("{name}.fmsrst"));
    let (straight_file, part_file, resumed_file) =
        (saved("straight"), saved("part"), saved("resumed"));
    std::thread::scope(|s| {
        s.spawn(|| {
            run(
                &dir.join("straight"),
                &[
                    &schedule[..],
                    &["--years", "4", "--save-restart", arg(&straight_file)],
                ]
                .concat(),
            );
        });
        let stderr = run(
            &dir.join("part"),
            &[
                &schedule[..],
                &["--years", "1.5", "--save-restart", arg(&part_file)],
            ]
            .concat(),
        );
        // Too short for its spin-up, but saved to continue.
        assert!(
            stderr.contains("warning: --years must be at least 4"),
            "{stderr}"
        );
    });
    let stderr = run(
        &dir.join("resumed"),
        &[
            &schedule[..],
            &[
                "--restart",
                arg(&part_file),
                "--years",
                "2.5",
                "--save-restart",
                arg(&resumed_file),
            ],
        ]
        .concat(),
    );
    assert!(stderr.contains("an exact continuation"), "{stderr}");
    assert!(!stderr.contains("warning"), "{stderr}");
    assert_same(&straight_file, &resumed_file);
}

#[test]
fn restarts_that_do_not_fit_are_refused() {
    let dir = TempDir::new("refused");
    let file = dir.join("small.fmsrst");
    run(
        &dir.join("small"),
        &[
            "--ocean-jumps",
            "0",
            "--years",
            "1.05",
            "--save-restart",
            arg(&file),
        ],
    );
    let not_restart = dir.join("notes.txt");
    std::fs::write(&not_restart, "not a restart").unwrap();
    let missing = dir.join("missing.fmsrst");
    for (path, frequency, reason) in [
        (&file, "5", "saved on a different grid"),
        (&not_restart, "4", "not a restart file"),
        (&missing, "4", "restart file:"),
    ] {
        let out = dir.join("refused");
        let args = ["--ocean-jumps", "0", "--years", "1", "--restart", arg(path)];
        let (ok, stderr) = climate(frequency, &out, &args);
        assert!(!ok, "{}: {stderr}", path.display());
        let expected = format!("error: cannot restart from {}: ", path.display());
        assert!(
            stderr.contains(&expected) && stderr.contains(reason),
            "{}: {stderr}",
            path.display()
        );
        assert!(!out.join("climate").exists(), "a refused run wrote output");
        eprintln!("{}", stderr.lines().last().unwrap());
    }
}
