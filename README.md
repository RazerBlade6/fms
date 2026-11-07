# Fantasy Map Simulator

A desktop app for generating, editing and simulating fantasy planets. You sculpt
terrain and water on a spherical world. A simplified climate model then produces
winds, ocean currents, weather and climate zones, and runs them forward in time.

Written in Rust (egui + wgpu). **Status:** early development (milestone 1).

```sh
cargo test --workspace      # run all tests
cargo run -p fms-cli        # headless runner (`fms`)
cargo run -p fms-app        # GUI (`fantasy-map`)
```

## Documentation

- [Design plan](docs/PLAN.md): architecture, key decisions, subsystems, equations
- [Roadmap](docs/ROADMAP.md): milestones and progress
