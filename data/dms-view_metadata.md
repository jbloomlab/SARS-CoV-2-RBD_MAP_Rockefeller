## Mutational escape from the polyclonal antibody response to SARS-CoV-2 infection is largely shaped by a single class of antibodies

For experimental background, see our paper **here (update with link)**.

### What data are shown here?
We are showing mutations to the SARS-CoV-2 RBD that escape binding by monoclonal antibodies or polyclonal convalescent human plasma, measured using mutational antigenic profiling. Raw data are available [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_MAP_Rockefeller/blob/main/results/supp_data/NY_mabs_sera_raw_data.csv).
The drop-down menus can be used to select the escape-mutation maps for each antibody or plasma.

When you click on sites, they will be highlighted on the protein structure of the ACE2-bound RBD ([PDB 6M0J](https://www.rcsb.org/structure/6M0J)) or the Fab-bound spike trimer or Fab-bound isolated RBD.

At the site level you can visualize one of two quantities:

 - *total escape* is the sum of the escape from all mutations at a site.

 - *max escape* is the magnitude of the largest-effect escape mutation at each site.

At the mutation level, the height of each letter is proportional to the extent to which that amino-acid mutation escapes antibody binding.
You can color the logo plot letters in four ways:

 - *escape color ACE2 bind* means color letters according to how that mutation affects ACE2 binding as measured in our prior deep mutational scanning ([Starr et al. (2020)](https://doi.org/10.1016/j.cell.2020.08.012), with yellow meaning highly deleterious, and brown meaning neutral or beneficial for ACE2 binding.

 - *escape color RBD expr* means color letters according to how that mutation affects RBD expression as measured in [Starr et al. (2020)](https://doi.org/10.1016/j.cell.2020.08.012).

 - *escape color gray* means color all letters gray.

 - *escape color func group* means color letters by functional group.
