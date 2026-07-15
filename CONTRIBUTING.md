# Contributing

Thanks for your interest in improving this workflow! This guide covers how to
report problems, set up a development environment, run the tests, and propose
changes.

By participating, you agree to abide by the [Code of Conduct](CODE_OF_CONDUCT.md).

## Reporting issues

Open a [GitHub issue](https://github.com/gynecoloji/snakemake_ATACseq/issues)
with:

- what you ran (the exact `snakemake` command / target),
- what you expected vs. what happened,
- the relevant log(s) from `logs/`, and your OS + Snakemake version.

## Development setup

```bash
git clone https://github.com/gynecoloji/snakemake_ATACseq.git
cd snakemake_ATACseq

# Driver environment (per-rule tool envs are created on first --use-conda run)
mamba create -n atacseq -c conda-forge -c bioconda snakemake pandas pytest
conda activate atacseq
```

Reference genomes/annotations are not tracked; download them into `ref/` as
described in [`config/README.md`](config/README.md).

## Running the workflow

```bash
snakemake -s workflow/Snakefile -n                       # dry run (validates config + DAG)
snakemake -s workflow/Snakefile --use-conda --cores 20   # everything (primary → QC)
snakemake -s workflow/Snakefile --use-conda --cores 20 atacseq_all   # primary only
```

Configuration is validated against [`workflow/schemas/config.schema.yaml`](workflow/schemas/config.schema.yaml)
on every run — that schema is the single source of truth for parameters.

## Tests

Unit tests (Python scripts) need only `numpy` + `pandas` and run fast:

```bash
python -m pytest tests/ -q
```

CI (`.github/workflows/ci.yml`) additionally builds the full DAG via a dry run
over stubbed inputs. Please make sure both pass before opening a PR.

## Commit messages & releases

This repo uses **[Conventional Commits](https://www.conventionalcommits.org)**
and [release-please](https://github.com/googleapis/release-please) to automate
versioning, the changelog, and releases. Prefix your commits:

| Prefix | Effect |
|---|---|
| `feat: …` | new feature → minor version bump, listed under *Added* |
| `fix: …` | bug fix → patch bump, under *Fixed* |
| `feat!: …` or a `BREAKING CHANGE:` footer | major bump |
| `docs:` / `refactor:` / `test:` / `chore:` | no release on their own |

On merge to `main`, release-please opens/updates a "release PR"; merging that PR
tags the version, publishes a GitHub Release, updates `CHANGELOG.md` and
`CITATION.cff`, and (via the Zenodo integration) mints a new DOI. You do **not**
edit the changelog or version numbers by hand.

## Pull requests

1. Branch from `main`, make your change, and add/adjust tests where relevant.
2. Run `pytest tests/ -q` and `snakemake -s workflow/Snakefile -n`.
3. Open a PR with a clear, Conventional-Commit-style title and description.
