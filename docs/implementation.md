# Implementation of Polaris

The software is written as a Python package, following current conventions (using a `src/polaris` layout and a `pyproject.toml` file.

## Runnable modules

Two of the modules in the package, the `polaris.pipeline` and `polaris.simulation` modules, are available as runnable Python modules, effectively providing two scripts: one for running the pipeline, and one for creating simulated (test) data. Technically, these can be easily turned into standalone scripts, but this has not been implemented to avoid cluttering the user's environment.

To run the modules, use the following

```
python -m polaris.pipeline
```

to run the actual pipeline, or

```
python -m polaris.simulation
```

to create simulated data.

Note that both modules have one required argument, a configuration file.

If you prefer standalone scripts, then start from the downloaded source code, uncomment the two lines below the `[project.script]` header in the `pyproject.toml` file (renaming the scripts as you prefer), and re-install the package with `pip install .` (or `pip install -e '.[dev]'` if you are actively developing the package).


## Configuration

Options are provided through a configuration file. This avoids dozens (eventually) of options, but can hide some configuration settings, as they are less explicit.

The configuration file format is TOML (Tom's Obvious Minimal Language), which provides a relatively easy to read configuration style, while keeping some extra configuration options (such as lists or dicts values) and checks (since it distinguishes integers, floating point values, date-time values and booleans).

The configuration files are currently provided as files in the root directory of the package, as `polaris.cfg` and `sim.cfg`.

The default pipelines are annotated with comments above a setting for a brief explanation of that setting. For more details, see the separate [pipeline](pipeline.md) and [simulation](simulation.md) documents.

Note that the two configuration files can be combined into one file, since the simulation part is contained inside its own section and subsection (or table and subtables in TOML language


#### Todo

The default configuration for each module can be printed to the standard output using the `--print-config` option. This output can then be redirected to a file that serves as a starting point for the configuration input when running the module, e.g. 

```
python -m polaris.pipeline --print-config > polaris.cfg
```
