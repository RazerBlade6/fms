# Fantasy Map Simulator

A desktop app for generating, editing and simulating fantasy planets. You sculpt
terrain and water on a spherical world. A simplified climate model then produces
winds, ocean currents, weather and climate zones, and runs them forward in time.

Written in Rust (egui + wgpu). **Status:** early development (M1 and M1.5 complete; M2 in progress).

```sh
cargo test --workspace                        # run all tests
cargo run --release -p fms-cli -- dump --all  # render every field to output/
cargo run -p fms-app                          # GUI (`fantasy-map`, placeholder)
```

## Documentation

- [Design plan](docs/PLAN.md): architecture, key decisions, subsystems, equations
- [Roadmap](docs/ROADMAP.md): milestones and progress
- [CLI](docs/CLI.md): the headless `fms` runner and its image output
- [Benchmarks](docs/BENCHMARKS.md): how to run them, and baseline numbers
- [Development](docs/DEVELOPMENT.md): setup, git workflow, determinism rules, conventions
- [Observations](docs/OBSERVATIONS.md): known limitations and ideas parked for later
