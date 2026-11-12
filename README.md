# Fantasy Map Simulator (fms)

A project centred around **generating, editing, simulating, and playing**
with fictional or theoretical maps on theoretical (roughly earth-like) planets.

Planets are assumed to be earth-like (so, earth-like atmospheric and lithic composition), and extreme
planetary states (tidally-locked worlds, highly eccentric orbits, binary or trinary star systems, etc.)
are not supported.

Most other things are customizable.

The project is divded into 4 modules:

- `fms-editor`: A Rust GUI App for editing and designing custom planets and terrain. Not built.
- `fms-climate`: The atmospheric GCM + Vegetation + Hydrological model, generating biomes, weather, water,
  currently in progress.
- `fms-geo`: A Geological/Tectonic plate simulator for generating physical and relief features,
  continental drift, carbonate-silicate cycle.
- `fms-civ`: A civilisation simulator, with factions, culture, religion, and probably more (Also not designed yet).

Written in Rust (egui + wgpu).

**Status:** Currently building `fms-climate`; M2 (the headless GCM) is in progress; fully built, requires tuning now.

This code is licensed under the PolyForm Non-Commercial License; Free for personal, hobbyist or academic use.
Commercial Use (for example, in a video game or a paid service) requires permission.

Since currently only climate works, this command will run the fms-climate module

```sh
cargo test --workspace                             # run all tests
cargo run --release -p fms-climate -- dump --all   # render every field to output/
cargo run -p fms-viewer                            # viewer (placeholder)
```

## Documentation

- [Plan](docs/PLAN.md): the `fms` family
- [Climate plan](docs/climate/PLAN.md): goals and key decisions of `fms-climate`
- [Roadmap](docs/climate/ROADMAP.md): milestones and progress
- [Design](docs/climate/DESIGN.md): how each subsystem works
- [Equations](docs/climate/EQUATIONS.tex): the model's equations in LaTeX, for review
- [`.fmsmap`](docs/FMSMAP.md): the file every module exchanges
- [CLI](docs/climate/CLI.md): the headless `fms-climate` runner and its image output
- [Benchmarks](docs/climate/BENCHMARKS.md): how to run them, and baseline numbers
- [Observations](docs/climate/OBSERVATIONS.md): known limitations and ideas parked for later
- [References](REFERENCES.md): the papers and data the design cites, with DOIs
