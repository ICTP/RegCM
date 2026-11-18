# rcmeval: RegCM evaluation against observations

`rcmeval` compares one or more RegCM simulations with observations and
reanalyses, and writes bias maps, region tables, region profiles on pressure
levels and a CSV file with every statistic.

    python -m rcmeval -c examples/sam_month.yaml

Requirements: numpy, scipy, xarray, dask, netCDF4, matplotlib, PyYAML;
cartopy for coastlines on the maps (optional). No ESMF is needed.

## How the comparison is made

* **Common grid.** The model and a reference are compared on the coarser of
  the two grids. The finer field is aggregated onto it: each fine cell goes
  to the coarse cell with the nearest centre, and a coarse value is the area
  weighted mean of its fine cells (conservative for fine to coarse within the
  accuracy of the cell assignment). Coarse cells less than half covered are
  left out. The weights are computed once and cached.
* **Order of operations.** Monthly means are computed once on the native
  grids, seasonal means from the monthly means (weighted by the month
  lengths), and only these are remapped. The 99th percentile is computed from
  the daily fields after remapping them to the common grid, so that model and
  reference extremes refer to the same scale.
* **Boundary.** `boundary` cells of the lateral relaxation zone are left out.
* **Surface.** Regions can be restricted to land or ocean (model land
  fraction on the common grid), and land-only (or ocean-only) datasets are
  compared only where they have data.
* **Statistics.** For every variable, season, region, run and reference:
  area weighted means, bias, RMSE and pattern correlation.
* **Cache.** Every intermediate result is stored under a key made of the
  source files (path, size, modification time), the period, the variable,
  the grids and the code of the package. A new model output, a changed
  observation file or a new version of `rcmeval` gives new keys: the cache
  never needs to be cleaned by hand (`--clean-cache` frees the space).

## Configuration

See `examples/sam_month.yaml`. The runs are RegCM output directories
(`format: regcm`, files `DOMAIN_TYPE.YYYYMMDDHH.nc`) or CORDEX archives
written by pycordexer (`format: cordex`). Each variable lists its references
by catalog name; a reference may use other years than the simulation, e.g. a
climatology:

    hfls: {refs: [ERA5, {name: OAFLUX, years: 2006-2015}]}

Output in `output`:

* `bias_VAR_RUN.png`: bias maps, seasons by row, references by column;
* `p99_VAR_RUN.png`: bias of the 99th percentile of daily values;
* `profile_VAR_REGION_SEASON.png`: region mean profiles;
* `cycle_VAR_RUN_REGION.png`: annual cycles (three months or more);
* `report.md`: tables of the region means and of bias, RMSE, correlation;
* `stats.csv`: all the statistics, for further processing.

## Catalogs of reference data

A catalog (`catalog/ictp.yaml`, `catalog/leonardo.yaml`) tells where the
datasets are on a machine and how to read them (file template, variable
name, conversion to the units of `rcmeval/units.py`). The builder writes one
for your machine, and tells you what is missing and where to download it:

    python -m rcmeval.builder --scan /path/to/OBS --scan /path/to/ERA5 -o mysite.yaml

or, for a remote machine, from a list of its files:

    ssh remote 'cd /path/to/OBS; find . -name "*.nc*"' > files.txt
    python -m rcmeval.builder --list files.txt --root /path/to/OBS -o remote.yaml

The datasets known to the builder, their file layouts and download sites are
in `rcmeval/known.yaml`; directories it does not recognize are listed, to be
added by hand to the site catalog or to `known.yaml`.

## Adding a variable or a dataset

* a new variable needs its target unit in `rcmeval/units.py`;
* a new dataset needs an entry in the site catalog (and preferably in
  `rcmeval/known.yaml`, with its download site, for the builder).
