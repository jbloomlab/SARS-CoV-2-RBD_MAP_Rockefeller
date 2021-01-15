# Custom analyses

This directory contains custom analyses for specific sets of antibodies or sera.

## Analyses in this directory:
### For Bjorkman / Barnes antibodies:
#### Write pymol commmands to make PSEs
The notebook [write_pymol_commands.ipynb](./write_pymol_commands.ipynb) takes as input the pdb files output by `output_pdbs.ipynb`, the escape sites called by `call_strong_escape_sites.ipynb`, and the contact sites called by `annotate_structural_contacts.Rmd`.

It writes text files that contain commands that can be pasted into the PyMol command line to create a PSE that can be used for further structure-gazing and manual analysis.

This notebook might work for other antibodies, but I haven't tested it.

#### Write pymol commmands to make 6M0J-based PSE
The notebook [pymol_commands_6m0j.ipynb](./pymol_commands_6m0j.ipynb) takes as input the pdb files output by `output_pdbs.ipynb`, the escape sites called by `call_strong_escape_sites.ipynb`, and the contact sites called by `annotate_structural_contacts.Rmd`.

It writes text files that contain commands that can be pasted into the PyMol command line to create a PSE that can be used for further structure-gazing and manual analysis, specifically with the 6M0J structure.

This notebook might work for other antibodies, but I haven't tested it.

#### Bjorkman custom analyses
The notebook [bjorkman_analyses.ipynb](bjorkman_analyses.ipynb) does some custom analyses for the Bjorkman antibodies and NY sera. 

### `custom-plots_xxx` scripts
These scripts make customized plots for different antibody sets that vary from the primary computational pipeline. They are currently just ran manually in R, but when each `xxx` antibody set splits off to its own repo, these can be incorporated into the `Snakemake` pipeline.