# Input data

This directory contains input data used for the analysis and configuration for some of the analyses:

 - [barcode_runs.csv](barcode_runs.csv): Illumina sequencing runs for barcode counting.
   The *sample* column should be the hyphen separated concatenation of *experiment*, *antibody*, *concentration*, and *selection*.
   The **combination** of the *library* and *sample* columns should be unique.
   The *frac_escape* column gives the overall fraction of the library that escaped antibody.
   The *R1* FASTQ files should be semicolon (`; `) separated list of glob patterns.

 - [wildtype_sequence.fasta](wildtype_sequence.fasta): wildtype sequence of mutagenized gene.

 - [./pdbs/](pdbs): files downloaded from the [Protein Data Bank](https://www.rcsb.org/) for structural analyses. The 7c01.pdb structure was modified to give the 7c01_single.pdb: the original structure has two Fab:RBD complexes within the asymmetric unit, so we removed chains B, C, and D to leave just one complex for visualization.

 - [escape_profiles_config.yaml](escape_profiles_config.yaml): Information on how to plot the escape profiles; manually edit this to alter their plotting.

 - [mds_config.yaml](mds_config.yaml): Information on how to do the multi-dimensional scaling.

 - [output_pdbs_config.yaml](output_pdbs_config.yaml): Information on how to output structural mappings of escape.

 - [structural_annotation_config.yaml](structural_annotation_config.yaml): Information on how to output structural annotations of contact residues in various RBD:ligand complexes.

 - [site_color_schemes.csv](site_color_schemes.csv): Schemes for how to color sites (can be used in escape profiles). Here are details on these schemes.

   - The `subdomain` scheme colors sites orange if they directly contact ACE2 (within 4 angstroms in PDB 6m0j) in the SARS-CoV-2 structure (residues 417, 446, 449, 453, 455, 456, 475, 486, 487, 489, 493, 496, 498, 500, 501, 502, 505), blue if they are in the receptor binding motif (RBM, residue 437 to 508, inclusive), and green if they are in the core RBD domain (all other sites). These definitions match those used in [Starr et al (2020)](https://www.cell.com/cell/fulltext/S0092-8674(20)31003-5).

 - [./RBD_sites.csv](RBD_sites.csv): Useful numeric and functional annotations of RBD sites for referencing in various analyses.

 - [./RBDs_aligned.fasta](RBDs_aligned.fasta): alignment of sarbecovirus RBDs

 - [escape_selection_results.yaml](escape_selection_results.yaml): results of viral escape-mutant selections, used to make plots about these.

 - [./Weisblum_SinoBiological_Spike.gb](Weisblum_SinoBiological_Spike.gb) is the spike sequence used in the [Weisblum et al. eLife (2020)](https://elifesciences.org/articles/61312) _in vitro_ VSV-Spike selection experiments. From the methods:

   * "The generation of infectious rVSV/SARS-CoV-2/GFP chimeric viruses
     stocks has been previously described (Schmidt et al., 2020). Two plaque purified
     variants designated rVSV/SARS-CoV-2/GFP1D7 and rVSV/SARS-CoV-2/
     GFP2E1 that encode F157S/R685M (1D7) and D215G/R683G (2E1)
     substitutions were used in these studies."

   * "For selection of viruses resistant to plasma or monoclonal
     antibodies, rVSV/SARS-CoV-2/GFP1D7 and rVSV/SARS-CoV-2/GFP2E1 
     populations containing 1e6 infectious particles were incubated with
     dilutions of monoclonal antibodies (10μg/ml, 5μg/ml) or COVID19
     plasma (1:50, 1:250, 1:500) for 1h at 37°C. "

   * [Schmidt, et al. (2020)](https://rupress.org/jem/article/217/11/e20201181/151961/Measuring-SARS-CoV-2-neutralizing-antibody) states:
     "To construct a replication competent rVSV/SARS-CoV-2 chimeric virus clone, a codon optimized cDNA sequence encoding the SARS-CoV-2 spike protein (SinoBiological) but lacking the C-terminal 18 codons was inserted, using Gibson cloning, into a recombinant VSV background that contains GFP immediately upstream of the L (polymerase) following a strategy we previously described for the exchange of VSV-G with HIV-1 Env proteins (Liberatore et al., 2019)"

   * Thus, we use the SinoBiological sequence from [here](https://www.sinobiological.com/cdna-clone/2019-ncov-cov-spike-vg40589-ut), called "SARS-CoV-2 (2019-nCoV) Spike ORF mammalian expression plasmid (Codon Optimized)."

  - [./pse_config.yaml](pse_config.yaml): config file for writing pymol commands to create a PSE session for certain selections.

  - [./pse_config_6m0j.yaml](pse_config_6m0j.yaml): config file for writing pymol commands to create a PSE session aligned to 6M0J (not ab-bound structures) for certain selections.

  - [./mds_color_schemes.csv](mds_color_schemes.csv): files that specifies point colors for MDS plots, rather than coloring pie charts by site.
