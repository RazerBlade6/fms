# fms-climate

A climate-centred simulator of theoretical planets, and part of the `fms` family of
modules that exchange `.fmsmap` files. Give it a planet and its terrain. A simplified
climate model then produces winds, ocean currents, weather, climate zones and
vegetation, and runs them forward in time.

Written in Rust (egui + wgpu). **Status:** early development (M1 and M1.5 complete; M2 in progress).

```sh
cargo test --workspace                             # run all tests
cargo run --release -p fms-climate -- dump --all   # render every field to output/
cargo run -p fms-viewer                            # viewer (placeholder)
```

## Documentation

- [Plan](docs/PLAN.md): the `fms` family, goals and key decisions
- [Roadmap](docs/ROADMAP.md): milestones and progress
- [Design](docs/DESIGN.md): how each subsystem works
- [Equations](docs/EQUATIONS.tex): the model's equations in LaTeX, for review
- [`.fmsmap`](docs/FMSMAP.md): the file every module exchanges
- [CLI](docs/CLI.md): the headless `fms-climate` runner and its image output
- [Benchmarks](docs/BENCHMARKS.md): how to run them, and baseline numbers
- [Development](docs/DEVELOPMENT.md): setup, git workflow, determinism rules, conventions
- [Observations](docs/OBSERVATIONS.md): known limitations and ideas parked for later
- [References](REFERENCES.md): the papers and data the design cites, with DOIs
