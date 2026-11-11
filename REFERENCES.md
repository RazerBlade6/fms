# References

Papers, books and datasets cited in the code and docs, with what we use them
for. Keep this list in step with the citations: add an entry when you cite a
new source. Each DOI was checked against a publisher, repository or
bibliographic record (November 2026); where no DOI exists or none could be
confirmed, the entry says so.

## Grid and numerics

- **Thuburn, J., Ringler, T. D., Skamarock, W. C. & Klemp, J. B. (2009).**
  Numerical representation of geostrophic modes on arbitrarily structured
  C-grids. *J. Comput. Phys.* 228, 8321–8335.
  [doi:10.1016/j.jcp.2009.08.006](https://doi.org/10.1016/j.jcp.2009.08.006)
  — TRiSK weights for the Coriolis term (`crates/grid`).
- **Ringler, T. D., Thuburn, J., Klemp, J. B. & Skamarock, W. C. (2010).**
  A unified approach to energy conservation and potential vorticity dynamics
  for arbitrarily-structured C-grids. *J. Comput. Phys.* 229, 3065–3090.
  [doi:10.1016/j.jcp.2009.12.007](https://doi.org/10.1016/j.jcp.2009.12.007)
  — the energy-conserving TRiSK operators and shallow-water scheme.
- **Perot, B. (2000).** Conservation properties of unstructured staggered mesh
  schemes. *J. Comput. Phys.* 159, 58–89.
  [doi:10.1006/jcph.2000.6424](https://doi.org/10.1006/jcph.2000.6424)
  — reconstructing vectors from edge normal components.
- **Peixoto, P. S. (2016).** Accuracy analysis of mimetic finite volume
  operators on geodesic grids and a consistent alternative. *J. Comput.
  Phys.* 310, 127–160.
  [doi:10.1016/j.jcp.2015.12.058](https://doi.org/10.1016/j.jcp.2015.12.058)
  — why the operators' errors don't converge at the pentagons
  ([OBSERVATIONS.md](docs/OBSERVATIONS.md) O-1).
- **Skilling, J. (2004).** Programming the Hilbert curve. *AIP Conf. Proc.*
  707, 381–387. (No DOI confirmed.) — the Hilbert-curve cell ordering.
- **Zalesak, S. T. (1979).** Fully multidimensional flux-corrected transport
  algorithms for fluids. *J. Comput. Phys.* 31, 335–362.
  [doi:10.1016/0021-9991(79)90051-2](https://doi.org/10.1016/0021-9991(79)90051-2)
  — the positive-definite humidity transport.
- **Robert, A. (1969).** The integration of a spectral model of the atmosphere
  by the implicit method. *Proc. WMO/IUGG Symposium on Numerical Weather
  Prediction*, Tokyo, VII-19–24. (No DOI.) — semi-implicit time stepping.
- **Simmons, A. J., Hoskins, B. J. & Burridge, D. M. (1978).** Stability of
  the semi-implicit method of time integration. *Mon. Wea. Rev.* 106,
  405–412.
  [doi:10.1175/1520-0493(1978)106<0405:SOTSIM>2.0.CO;2](https://doi.org/10.1175/1520-0493(1978)106<0405:SOTSIM>2.0.CO;2)
  (DOI from citing papers' reference lists.) — why the semi-implicit
  reference state is isothermal and warm.
- **Williams, P. D. (2009).** A proposed modification to the Robert–Asselin
  time filter. *Mon. Wea. Rev.* 137, 2538–2546.
  [doi:10.1175/2009MWR2724.1](https://doi.org/10.1175/2009MWR2724.1)
  — the RAW filter on the leapfrog steps.

## Dynamical core and benchmarks

- **Phillips, N. A. (1954).** Energy transformations and meridional
  circulations associated with simple baroclinic waves in a two-level,
  quasi-geostrophic model. *Tellus* 6, 273–286.
  [doi:10.1111/j.2153-3490.1954.tb01123.x](https://doi.org/10.1111/j.2153-3490.1954.tb01123.x)
  — the two-level model and its baroclinic instability.
- **Arakawa, A. & Suarez, M. J. (1983).** Vertical differencing of the
  primitive equations in sigma coordinates. *Mon. Wea. Rev.* 111, 34–45.
  [doi:10.1175/1520-0493(1983)111<0034:VDOTPE>2.0.CO;2](https://doi.org/10.1175/1520-0493(1983)111<0034:VDOTPE>2.0.CO;2)
  — the energy-conserving vertical discretization, for any number of layers.
- **Holton, J. R.** *An Introduction to Dynamic Meteorology.* Academic Press.
  (Textbook.) — two-level quasi-geostrophic theory (§8.2).
- **Williamson, D. L., Drake, J. B., Hack, J. J., Jakob, R. &
  Swarztrauber, P. N. (1992).** A standard test set for numerical
  approximations to the shallow water equations in spherical geometry.
  *J. Comput. Phys.* 102, 211–224.
  [doi:10.1016/S0021-9991(05)80016-6](https://doi.org/10.1016/S0021-9991(05)80016-6)
  — shallow-water test cases 2, 5 and 6.
- **Held, I. M. & Suarez, M. J. (1994).** A proposal for the intercomparison
  of the dynamical cores of atmospheric general circulation models. *Bull.
  Amer. Meteor. Soc.* 75, 1825–1830.
  [doi:10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2](https://doi.org/10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2)
  — the dry dynamical-core benchmark.

## Radiation

- **Manabe, S. & Wetherald, R. T. (1967).** Thermal equilibrium of the
  atmosphere with a given distribution of relative humidity. *J. Atmos.
  Sci.* 24, 241–259.
  [doi:10.1175/1520-0469(1967)024<0241:TEOTAW>2.0.CO;2](https://doi.org/10.1175/1520-0469(1967)024<0241:TEOTAW>2.0.CO;2)
  — the relative-humidity profile shape.
- **Myhre, G., Highwood, E. J., Shine, K. P. & Stordal, F. (1998).** New
  estimates of radiative forcing due to well mixed greenhouse gases.
  *Geophys. Res. Lett.* 25, 2715–2718.
  [doi:10.1029/98GL01908](https://doi.org/10.1029/98GL01908)
  — the CO₂ forcing the radiation is calibrated to.
- **Wilson, D. J. & Gea-Banacloche, J. (2012).** Simple model to estimate the
  contribution of atmospheric CO₂ to the Earth's greenhouse effect. *Am. J.
  Phys.* 80, 306.
  [doi:10.1119/1.3681188](https://doi.org/10.1119/1.3681188)
  — the CO₂ band-wing model.
- **Jeevanjee, N., Seeley, J. T., Paynter, D. & Fueglistaler, S. (2021).** An
  analytical model for spatially varying clear-sky CO₂ forcing. *J.
  Climate* 34, 9463–9480.
  [doi:10.1175/JCLI-D-19-0756.1](https://doi.org/10.1175/JCLI-D-19-0756.1)
  — logarithmic CO₂ forcing from the band edges.
- **Alduchov, O. A. & Eskridge, R. E. (1996).** Improved Magnus form
  approximation of saturation vapor pressure. *J. Appl. Meteor.* 35,
  601–609.
  [doi:10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2](https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2)
  — the saturation vapour pressure over water and ice.

## Moisture, convection and clouds

- **Manabe, S., Smagorinsky, J. & Strickler, R. F. (1965).** Simulated
  climatology of a general circulation model with a hydrologic cycle. *Mon.
  Wea. Rev.* 93, 769–798.
  [doi:10.1175/1520-0493(1965)093<0769:SCOAGC>2.3.CO;2](https://doi.org/10.1175/1520-0493(1965)093<0769:SCOAGC>2.3.CO;2)
  — moist convective adjustment.
- **Betts, A. K. (1986).** A new convective adjustment scheme. Part I:
  Observational and theoretical basis. *Quart. J. Roy. Meteor. Soc.* 112,
  677–691.
  [doi:10.1002/qj.49711247307](https://doi.org/10.1002/qj.49711247307)
  — Betts–Miller convection.
- **Frierson, D. M. W. (2007).** The dynamics of idealized convection schemes
  and their effect on the zonally averaged tropical circulation. *J. Atmos.
  Sci.* 64, 1959–1976.
  [doi:10.1175/JAS3935.1](https://doi.org/10.1175/JAS3935.1)
  — the simplified Betts–Miller scheme (deep and shallow, energy
  conserving).
- **Molteni, F. (2003).** Atmospheric simulations using a GCM with simplified
  physical parametrizations. I: Model climatology and variability in
  multi-decadal experiments. *Clim. Dyn.* 20, 175–191.
  [doi:10.1007/s00382-002-0268-2](https://doi.org/10.1007/s00382-002-0268-2)
  — SPEEDY: clouds from relative humidity and convective rain.
- **Vavrus, S. & Waliser, D. (2008).** An improved parameterization for
  simulating Arctic cloud amount in the CCSM3 climate model. *J. Climate.*
  [doi:10.1175/2008JCLI2299.1](https://doi.org/10.1175/2008JCLI2299.1)
  — "freeze-dry": low cloud thins in air below 3 g/kg of water vapour.

## Surface and boundary layer

- **Lott, F. & Miller, M. J. (1997).** A new subgrid-scale orographic drag
  parametrization: Its formulation and testing. *Quart. J. Roy. Meteor.
  Soc.* 123, 101–127.
  [doi:10.1002/qj.49712353704](https://doi.org/10.1002/qj.49712353704)
  — planned sub-grid orographic drag (M2 step 10).

## Ocean and sea ice

- **Stommel, H. (1948).** The westward intensification of wind-driven ocean
  currents. *Trans. Amer. Geophys. Union* 29, 202–206.
  [doi:10.1029/TR029i002p00202](https://doi.org/10.1029/TR029i002p00202)
  — the wind-driven gyre and western boundary currents (M2 step 12).
- **Munk, W. H. (1950).** On the wind-driven ocean circulation. *J. Meteor.*
  7, 80–93.
  [doi:10.1175/1520-0469(1950)007<0080:OTWDOC>2.0.CO;2](https://doi.org/10.1175/1520-0469(1950)007<0080:OTWDOC>2.0.CO;2)
  — the gyre with lateral friction (M2 step 12).
- **Gent, P. R. & McWilliams, J. C. (1990).** Isopycnal mixing in ocean
  circulation models. *J. Phys. Oceanogr.* 20, 150–155.
  [doi:10.1175/1520-0485(1990)020<0150:IMIOCM>2.0.CO;2](https://doi.org/10.1175/1520-0485(1990)020<0150:IMIOCM>2.0.CO;2)
  — eddy heat transport in the slab ocean (M2 step 12).
- **Semtner, A. J. (1976).** A model for the thermodynamic growth of sea ice
  in numerical investigations of climate. *J. Phys. Oceanogr.* 6, 379–389.
  [doi:10.1175/1520-0485(1976)006<0379:AMFTTG>2.0.CO;2](https://doi.org/10.1175/1520-0485(1976)006<0379:AMFTTG>2.0.CO;2)
  — zero-layer thermodynamic sea ice (M2 step 11).
- **Codron, F. (2012).** Ekman heat transport for slab oceans. *Clim. Dyn.*
  38, 379–389.
  [doi:10.1007/s00382-011-1031-3](https://doi.org/10.1007/s00382-011-1031-3)
  — the two-layer slab ocean with Ekman transport (M2 steps 11–12).
- **Bhatnagar, S., Codron, F., Millour, E., Bolmont, E., Brunetti, M.,
  Kasparian, J., Turbet, M. & Chaverot, G. (2026).** A fast and physically
  grounded ocean model for GCMs: the Dynamical Slab Ocean Model of the
  Generic-PCM (rev. 3423). *Geosci. Model Dev.* 19, 3285–3316.
  [doi:10.5194/gmd-19-3285-2026](https://doi.org/10.5194/gmd-19-3285-2026)
  — the ocean design of M2 steps 11–12 (Sverdrup blend, Gent–McWilliams,
  sea ice albedo) and its aquaplanet benchmark (PLAN §4.3). Its code and run
  settings are archived at
  [doi:10.5281/zenodo.18771594](https://doi.org/10.5281/zenodo.18771594); we
  follow them where they differ from the paper's text.

## Planet, terrain, hydrology and climate classification

- **Berger, A. L. (1978).** Long-term variations of daily insolation and
  Quaternary climatic changes. *J. Atmos. Sci.* 35, 2362–2367.
  [doi:10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2](https://doi.org/10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2)
  — daily-mean insolation.
- **NOAA National Centers for Environmental Information (2022).** ETOPO 2022
  15 Arc-Second Global Relief Model.
  [doi:10.25921/fd45-gt74](https://doi.org/10.25921/fd45-gt74)
  — Earth's heightmap for the benchmarks.
- **Cortial, Y., Peytavie, A., Galin, É. & Guérin, É. (2019).** Procedural
  tectonic planets. *Comput. Graph. Forum* 38(2), 1–11.
  [doi:10.1111/cgf.13614](https://doi.org/10.1111/cgf.13614)
  — reading for `fms-geo`'s plate tectonics.
- **Barnes, R., Lehman, C. & Mulla, D. (2014).** Priority-flood: An optimal
  depression-filling and watershed-labeling algorithm for digital elevation
  models. *Comput. Geosci.* 62, 117–127.
  [doi:10.1016/j.cageo.2013.04.024](https://doi.org/10.1016/j.cageo.2013.04.024)
  — depression filling for rivers and lakes (M3).
- **Peel, M. C., Finlayson, B. L. & McMahon, T. A. (2007).** Updated world map
  of the Köppen–Geiger climate classification. *Hydrol. Earth Syst. Sci.*
  11, 1633–1644.
  [doi:10.5194/hess-11-1633-2007](https://doi.org/10.5194/hess-11-1633-2007)
  — the Köppen–Geiger criteria and the Earth map to score against (M3).

## Earth reference climatology

Kept in `reference/earth` (see its SOURCES.md) and compared against by
`scripts/compare_earth.py`.

- **Kalnay, E., Kanamitsu, M., Kistler, R., et al. (1996).** The NCEP/NCAR
  40-year reanalysis project. *Bull. Amer. Meteor. Soc.* 77, 437–471.
  [doi:10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2](https://doi.org/10.1175/1520-0477(1996)077%3C0437:TNYRP%3E2.0.CO;2)
  — monthly skin and 2 m temperature, 1991–2020 long-term means from NOAA PSL.
- **Adler, R. F., Sapiano, M. R. P., Huffman, G. J., et al. (2018).** The
  Global Precipitation Climatology Project (GPCP) monthly analysis (new
  version 2.3) and a review of 2017 global precipitation. *Atmosphere* 9, 138.
  [doi:10.3390/atmos9040138](https://doi.org/10.3390/atmos9040138)
  — monthly precipitation, 1991–2020 long-term means from NOAA PSL.
- **Fetterer, F., Knowles, K., Meier, W. & Savoie, M. (2002, updated
  daily).** Sea Ice Index, Version 4 (G02135). National Snow and Ice Data
  Center. [doi:10.7265/N5QJ7F7W](https://doi.org/10.7265/N5QJ7F7W) (the DOI
  the v4 user guide gives for the series; no version-4 DOI confirmed)
  — monthly sea-ice area and extent per hemisphere, and the Arctic pole-hole
  sizes (user guide, Table 8).
