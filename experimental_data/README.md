# Experimental data for manual analysis
This directory contains experimental data that are analyzed "manually" for validation experiments.
These are outside the make `Snakemake` workflow.

## Neutralization experiments from Weisblum, et al.

[Weisblum et al. (2020)](https://elifesciences.org/articles/61312) does really nice work testing the effects of many mutations on neutralization by C121, C135, and C144 monoclonal antibodies, as well as against COV-47, COV-72, COV-107, and COV-NY sera (the last of which we do not examine in our study).

The IC50 and NT50 values from this paper are contained in the ./data/ directory.

This directory contains:
- [./data/](./data/) - the IC50 and NT50 values from Weisblum et al.
- [./results/](./results/) - plots for potential figures
- [neutralization_plots.ipynb](neutralization_plots.ipynb) - replotting data Weisblum data as fold-change IC50 relative to wildtype spike, correlated with our escape measurements.
