# Output PDBs with escape scores as B factors
This Python Jupyter notebook outputs PDBs with the escape scores as B factors.

Though we will want more elaborate series of commands to codify our visualization of these RBD structures colored by escape, the series of commands below, when executed in a `PyMol` session with one of these PDBs open, will color the RBD surface according to escape scores.

For example, to normalize each structure colored by the max mut effect, we might want to have a white to red scale from 0 to 1:

     create RBD, chain E
     hide all; show cartoon, chain A; color gray20, chain A
     show surface, RBD; spectrum b, white red, RBD, minimum=0, maximum=1
     
For something like total escape, maybe we want each structure normalized to the maximum total escape in that structure, in which case we can just leave the maximum argument empty:

     create RBD, chain E
     hide all; show cartoon, chain A; color gray20, chain A
     show surface, RBD; spectrum b, white red, RBD, minimum=0
     
We write PDBs with B factors indicating the total site escape and maximum mutation escape at each site, and the same with these values normalized to the maximum for the full structure (the latter are easier to process in `Chimera`).

First, import Python modules:


```python
import collections
import copy
import os
import warnings

import Bio.PDB

import dms_variants.pdb_utils

from IPython.display import display, HTML

import pandas as pd

import yaml
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Read configuration for outputting PDBs:


```python
print(f"Reading PDB output configuration from {config['output_pdbs_config']}")
with open(config['output_pdbs_config']) as f:
    output_pdbs_config = yaml.safe_load(f)
```

    Reading PDB output configuration from data/output_pdbs_config.yaml


Make output directory:


```python
os.makedirs(config['pdb_outputs_dir'], exist_ok=True)
```

Read escape fractions and compute **total** and **maximum** escape at each site, and also the total and maximum escape at each site normalized to be between 0 and 1 for each selection:


```python
print(f"Reading escape fractions from {config['escape_fracs']}")

escape_fracs = (
    pd.read_csv(config['escape_fracs'])
    .query('library == "average"')
    .assign(site=lambda x: x['label_site'])
    .groupby(['selection', 'site'])
    .aggregate(total_escape=pd.NamedAgg(config['mut_metric'], 'sum'),
               max_escape=pd.NamedAgg(config['mut_metric'], 'max')
               )
    .reset_index()
    .assign(max_total_escape=lambda x: x.groupby('selection')['total_escape'].transform('max'),
            max_max_escape=lambda x: x.groupby('selection')['max_escape'].transform('max'),
            norm_total_escape=lambda x: x['total_escape'] / x['max_total_escape'],
            norm_max_escape=lambda x: x['max_escape'] / x['max_max_escape'])
    )

display(HTML(escape_fracs.head().to_html(index=False)))
```

    Reading escape fractions from results/escape_scores/escape_fracs.csv



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>selection</th>
      <th>site</th>
      <th>total_escape</th>
      <th>max_escape</th>
      <th>max_total_escape</th>
      <th>max_max_escape</th>
      <th>norm_total_escape</th>
      <th>norm_max_escape</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>0.049257</td>
      <td>0.007632</td>
      <td>2.99825</td>
      <td>0.2456</td>
      <td>0.016428</td>
      <td>0.031075</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>332</td>
      <td>0.109619</td>
      <td>0.015340</td>
      <td>2.99825</td>
      <td>0.2456</td>
      <td>0.036561</td>
      <td>0.062459</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>333</td>
      <td>0.051312</td>
      <td>0.008701</td>
      <td>2.99825</td>
      <td>0.2456</td>
      <td>0.017114</td>
      <td>0.035428</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>334</td>
      <td>0.153061</td>
      <td>0.030090</td>
      <td>2.99825</td>
      <td>0.2456</td>
      <td>0.051050</td>
      <td>0.122516</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>335</td>
      <td>0.115742</td>
      <td>0.021650</td>
      <td>2.99825</td>
      <td>0.2456</td>
      <td>0.038603</td>
      <td>0.088151</td>
    </tr>
  </tbody>
</table>


Now map the escape metrics to the B-factors.
For sites where no mutations have escape scores:
 - In the RBD chain(s) fill the B-factor for non-normalized scores to -1 to enable collapsing to zero or callout as a a separate class, depending how we choose to color sites for different visualizations. For normalized scores, fill to 0.
 - In other chains, always fill missing B factors to 0.  


```python
for name, specs in output_pdbs_config.items():
    print(f"\nMaking PDB mappings for {name} to {specs['pdbfile']}")
    assert os.path.isfile(specs['pdbfile'])
    
    # get escape fracs just for conditions of interest
    if isinstance(specs['conditions'], str) and specs['conditions'].upper() == 'ALL':
        conditions = escape_fracs['selection'].unique().tolist()
    else:
        assert isinstance(specs['conditions'], list)
        conditions = specs['conditions']
    print(f"Making mappings for {len(conditions)} conditions.")
    df = escape_fracs.query('selection in @conditions')
    
    # get chains
    assert isinstance(specs['chains'], list)
    print('Mapping to the following chains: ' + ', '.join(specs['chains']))
    df = pd.concat([df.assign(chain=chain) for chain in specs['chains']], ignore_index=True)
    
    # make mappings for each condition and metric
    for condition, df in df.groupby('selection'):
        print(f"  Writing B-factor re-assigned PDBs for {condition} to:")
    
        for metric in ['total_escape', 'max_escape', 'norm_total_escape', 'norm_max_escape']:
        
            # what do we assign to missing sites?
            missing_metric = collections.defaultdict(lambda: 0)  # non-RBD chains always fill to zero
            for chain in specs['chains']:
                if 'norm' in metric:
                    missing_metric[chain] = 0  # missing sites in RBD are 0 for normalized metric PDBs
                else:
                    missing_metric[chain] = -1  # missing sites in RBD are -1 for non-normalized metric PDBs
        
            fname = os.path.join(config['pdb_outputs_dir'], f"{condition}_{name}_{metric}.pdb")
            print(f"    {fname}")
            
            dms_variants.pdb_utils.reassign_b_factor(input_pdbfile=specs['pdbfile'],
                                                     output_pdbfile=fname,
                                                     df=df,
                                                     metric_col=metric,
                                                     missing_metric=missing_metric)
```

    
    Making PDB mappings for 6m0j to data/pdbs/6M0J.pdb
    Making mappings for 53 conditions.
    Mapping to the following chains: E
      Writing B-factor re-assigned PDBs for 12C_d152_80 to:
        results/pdb_outputs/12C_d152_80_6m0j_total_escape.pdb
        results/pdb_outputs/12C_d152_80_6m0j_max_escape.pdb
        results/pdb_outputs/12C_d152_80_6m0j_norm_total_escape.pdb
        results/pdb_outputs/12C_d152_80_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 12C_d61_160 to:
        results/pdb_outputs/12C_d61_160_6m0j_total_escape.pdb
        results/pdb_outputs/12C_d61_160_6m0j_max_escape.pdb
        results/pdb_outputs/12C_d61_160_6m0j_norm_total_escape.pdb
        results/pdb_outputs/12C_d61_160_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 13_d121_1250 to:
        results/pdb_outputs/13_d121_1250_6m0j_total_escape.pdb
        results/pdb_outputs/13_d121_1250_6m0j_max_escape.pdb
        results/pdb_outputs/13_d121_1250_6m0j_norm_total_escape.pdb
        results/pdb_outputs/13_d121_1250_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 13_d15_200 to:
        results/pdb_outputs/13_d15_200_6m0j_total_escape.pdb
        results/pdb_outputs/13_d15_200_6m0j_max_escape.pdb
        results/pdb_outputs/13_d15_200_6m0j_norm_total_escape.pdb
        results/pdb_outputs/13_d15_200_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 1C_d113_200 to:
        results/pdb_outputs/1C_d113_200_6m0j_total_escape.pdb
        results/pdb_outputs/1C_d113_200_6m0j_max_escape.pdb
        results/pdb_outputs/1C_d113_200_6m0j_norm_total_escape.pdb
        results/pdb_outputs/1C_d113_200_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 1C_d26_200 to:
        results/pdb_outputs/1C_d26_200_6m0j_total_escape.pdb
        results/pdb_outputs/1C_d26_200_6m0j_max_escape.pdb
        results/pdb_outputs/1C_d26_200_6m0j_norm_total_escape.pdb
        results/pdb_outputs/1C_d26_200_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 2196+2130_400 to:
        results/pdb_outputs/2196+2130_400_6m0j_total_escape.pdb
        results/pdb_outputs/2196+2130_400_6m0j_max_escape.pdb
        results/pdb_outputs/2196+2130_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/2196+2130_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 22C_d104_200 to:
        results/pdb_outputs/22C_d104_200_6m0j_total_escape.pdb
        results/pdb_outputs/22C_d104_200_6m0j_max_escape.pdb
        results/pdb_outputs/22C_d104_200_6m0j_norm_total_escape.pdb
        results/pdb_outputs/22C_d104_200_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 22C_d28_200 to:
        results/pdb_outputs/22C_d28_200_6m0j_total_escape.pdb
        results/pdb_outputs/22C_d28_200_6m0j_max_escape.pdb
        results/pdb_outputs/22C_d28_200_6m0j_norm_total_escape.pdb
        results/pdb_outputs/22C_d28_200_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 23C_d102_80 to:
        results/pdb_outputs/23C_d102_80_6m0j_total_escape.pdb
        results/pdb_outputs/23C_d102_80_6m0j_max_escape.pdb
        results/pdb_outputs/23C_d102_80_6m0j_norm_total_escape.pdb
        results/pdb_outputs/23C_d102_80_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 23C_d26_80 to:
        results/pdb_outputs/23C_d26_80_6m0j_total_escape.pdb
        results/pdb_outputs/23C_d26_80_6m0j_max_escape.pdb
        results/pdb_outputs/23C_d26_80_6m0j_norm_total_escape.pdb
        results/pdb_outputs/23C_d26_80_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 23_d120_500 to:
        results/pdb_outputs/23_d120_500_6m0j_total_escape.pdb
        results/pdb_outputs/23_d120_500_6m0j_max_escape.pdb
        results/pdb_outputs/23_d120_500_6m0j_norm_total_escape.pdb
        results/pdb_outputs/23_d120_500_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 23_d21_1250 to:
        results/pdb_outputs/23_d21_1250_6m0j_total_escape.pdb
        results/pdb_outputs/23_d21_1250_6m0j_max_escape.pdb
        results/pdb_outputs/23_d21_1250_6m0j_norm_total_escape.pdb
        results/pdb_outputs/23_d21_1250_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 23_d45_1250 to:
        results/pdb_outputs/23_d45_1250_6m0j_total_escape.pdb
        results/pdb_outputs/23_d45_1250_6m0j_max_escape.pdb
        results/pdb_outputs/23_d45_1250_6m0j_norm_total_escape.pdb
        results/pdb_outputs/23_d45_1250_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 24C_d104_200 to:
        results/pdb_outputs/24C_d104_200_6m0j_total_escape.pdb
        results/pdb_outputs/24C_d104_200_6m0j_max_escape.pdb
        results/pdb_outputs/24C_d104_200_6m0j_norm_total_escape.pdb
        results/pdb_outputs/24C_d104_200_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 24C_d32_200 to:
        results/pdb_outputs/24C_d32_200_6m0j_total_escape.pdb
        results/pdb_outputs/24C_d32_200_6m0j_max_escape.pdb
        results/pdb_outputs/24C_d32_200_6m0j_norm_total_escape.pdb
        results/pdb_outputs/24C_d32_200_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 25C_d115_80 to:
        results/pdb_outputs/25C_d115_80_6m0j_total_escape.pdb
        results/pdb_outputs/25C_d115_80_6m0j_max_escape.pdb
        results/pdb_outputs/25C_d115_80_6m0j_norm_total_escape.pdb
        results/pdb_outputs/25C_d115_80_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 25C_d48_200 to:
        results/pdb_outputs/25C_d48_200_6m0j_total_escape.pdb
        results/pdb_outputs/25C_d48_200_6m0j_max_escape.pdb
        results/pdb_outputs/25C_d48_200_6m0j_norm_total_escape.pdb
        results/pdb_outputs/25C_d48_200_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 25_d18_500 to:
        results/pdb_outputs/25_d18_500_6m0j_total_escape.pdb
        results/pdb_outputs/25_d18_500_6m0j_max_escape.pdb
        results/pdb_outputs/25_d18_500_6m0j_norm_total_escape.pdb
        results/pdb_outputs/25_d18_500_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 25_d94_200 to:
        results/pdb_outputs/25_d94_200_6m0j_total_escape.pdb
        results/pdb_outputs/25_d94_200_6m0j_max_escape.pdb
        results/pdb_outputs/25_d94_200_6m0j_norm_total_escape.pdb
        results/pdb_outputs/25_d94_200_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 6C_d33_500 to:
        results/pdb_outputs/6C_d33_500_6m0j_total_escape.pdb
        results/pdb_outputs/6C_d33_500_6m0j_max_escape.pdb
        results/pdb_outputs/6C_d33_500_6m0j_norm_total_escape.pdb
        results/pdb_outputs/6C_d33_500_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 6C_d76_500 to:
        results/pdb_outputs/6C_d76_500_6m0j_total_escape.pdb
        results/pdb_outputs/6C_d76_500_6m0j_max_escape.pdb
        results/pdb_outputs/6C_d76_500_6m0j_norm_total_escape.pdb
        results/pdb_outputs/6C_d76_500_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 7C_d103_200 to:
        results/pdb_outputs/7C_d103_200_6m0j_total_escape.pdb
        results/pdb_outputs/7C_d103_200_6m0j_max_escape.pdb
        results/pdb_outputs/7C_d103_200_6m0j_norm_total_escape.pdb
        results/pdb_outputs/7C_d103_200_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for 7C_d29_500 to:
        results/pdb_outputs/7C_d29_500_6m0j_total_escape.pdb
        results/pdb_outputs/7C_d29_500_6m0j_max_escape.pdb
        results/pdb_outputs/7C_d29_500_6m0j_norm_total_escape.pdb
        results/pdb_outputs/7C_d29_500_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for C002_400 to:
        results/pdb_outputs/C002_400_6m0j_total_escape.pdb
        results/pdb_outputs/C002_400_6m0j_max_escape.pdb
        results/pdb_outputs/C002_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/C002_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for C105_400 to:
        results/pdb_outputs/C105_400_6m0j_total_escape.pdb
        results/pdb_outputs/C105_400_6m0j_max_escape.pdb
        results/pdb_outputs/C105_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/C105_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for C110_400 to:
        results/pdb_outputs/C110_400_6m0j_total_escape.pdb
        results/pdb_outputs/C110_400_6m0j_max_escape.pdb
        results/pdb_outputs/C110_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/C110_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for C121_400 to:
        results/pdb_outputs/C121_400_6m0j_total_escape.pdb
        results/pdb_outputs/C121_400_6m0j_max_escape.pdb
        results/pdb_outputs/C121_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/C121_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for C135_400 to:
        results/pdb_outputs/C135_400_6m0j_total_escape.pdb
        results/pdb_outputs/C135_400_6m0j_max_escape.pdb
        results/pdb_outputs/C135_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/C135_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for C144_400 to:
        results/pdb_outputs/C144_400_6m0j_total_escape.pdb
        results/pdb_outputs/C144_400_6m0j_max_escape.pdb
        results/pdb_outputs/C144_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/C144_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for CB6_400 to:
        results/pdb_outputs/CB6_400_6m0j_total_escape.pdb
        results/pdb_outputs/CB6_400_6m0j_max_escape.pdb
        results/pdb_outputs/CB6_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/CB6_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV-021_500 to:
        results/pdb_outputs/COV-021_500_6m0j_total_escape.pdb
        results/pdb_outputs/COV-021_500_6m0j_max_escape.pdb
        results/pdb_outputs/COV-021_500_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV-021_500_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV-047_200 to:
        results/pdb_outputs/COV-047_200_6m0j_total_escape.pdb
        results/pdb_outputs/COV-047_200_6m0j_max_escape.pdb
        results/pdb_outputs/COV-047_200_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV-047_200_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV-057_50 to:
        results/pdb_outputs/COV-057_50_6m0j_total_escape.pdb
        results/pdb_outputs/COV-057_50_6m0j_max_escape.pdb
        results/pdb_outputs/COV-057_50_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV-057_50_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV-072_200 to:
        results/pdb_outputs/COV-072_200_6m0j_total_escape.pdb
        results/pdb_outputs/COV-072_200_6m0j_max_escape.pdb
        results/pdb_outputs/COV-072_200_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV-072_200_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV-107_80 to:
        results/pdb_outputs/COV-107_80_6m0j_total_escape.pdb
        results/pdb_outputs/COV-107_80_6m0j_max_escape.pdb
        results/pdb_outputs/COV-107_80_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV-107_80_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV2-2050_400 to:
        results/pdb_outputs/COV2-2050_400_6m0j_total_escape.pdb
        results/pdb_outputs/COV2-2050_400_6m0j_max_escape.pdb
        results/pdb_outputs/COV2-2050_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV2-2050_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV2-2082_400 to:
        results/pdb_outputs/COV2-2082_400_6m0j_total_escape.pdb
        results/pdb_outputs/COV2-2082_400_6m0j_max_escape.pdb
        results/pdb_outputs/COV2-2082_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV2-2082_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV2-2094_400 to:
        results/pdb_outputs/COV2-2094_400_6m0j_total_escape.pdb
        results/pdb_outputs/COV2-2094_400_6m0j_max_escape.pdb
        results/pdb_outputs/COV2-2094_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV2-2094_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV2-2096_400 to:
        results/pdb_outputs/COV2-2096_400_6m0j_total_escape.pdb
        results/pdb_outputs/COV2-2096_400_6m0j_max_escape.pdb
        results/pdb_outputs/COV2-2096_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV2-2096_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV2-2130_400 to:
        results/pdb_outputs/COV2-2130_400_6m0j_total_escape.pdb
        results/pdb_outputs/COV2-2130_400_6m0j_max_escape.pdb
        results/pdb_outputs/COV2-2130_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV2-2130_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV2-2165_400 to:
        results/pdb_outputs/COV2-2165_400_6m0j_total_escape.pdb
        results/pdb_outputs/COV2-2165_400_6m0j_max_escape.pdb
        results/pdb_outputs/COV2-2165_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV2-2165_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV2-2196_400 to:
        results/pdb_outputs/COV2-2196_400_6m0j_total_escape.pdb
        results/pdb_outputs/COV2-2196_400_6m0j_max_escape.pdb
        results/pdb_outputs/COV2-2196_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV2-2196_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV2-2479_400 to:
        results/pdb_outputs/COV2-2479_400_6m0j_total_escape.pdb
        results/pdb_outputs/COV2-2479_400_6m0j_max_escape.pdb
        results/pdb_outputs/COV2-2479_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV2-2479_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV2-2499_400 to:
        results/pdb_outputs/COV2-2499_400_6m0j_total_escape.pdb
        results/pdb_outputs/COV2-2499_400_6m0j_max_escape.pdb
        results/pdb_outputs/COV2-2499_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV2-2499_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV2-2677_400 to:
        results/pdb_outputs/COV2-2677_400_6m0j_total_escape.pdb
        results/pdb_outputs/COV2-2677_400_6m0j_max_escape.pdb
        results/pdb_outputs/COV2-2677_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV2-2677_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for COV2-2832_400 to:
        results/pdb_outputs/COV2-2832_400_6m0j_total_escape.pdb
        results/pdb_outputs/COV2-2832_400_6m0j_max_escape.pdb
        results/pdb_outputs/COV2-2832_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/COV2-2832_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for CR3022_400 to:
        results/pdb_outputs/CR3022_400_6m0j_total_escape.pdb
        results/pdb_outputs/CR3022_400_6m0j_max_escape.pdb
        results/pdb_outputs/CR3022_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/CR3022_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for LY-CoV555+LY-CoV016_400 to:
        results/pdb_outputs/LY-CoV555+LY-CoV016_400_6m0j_total_escape.pdb
        results/pdb_outputs/LY-CoV555+LY-CoV016_400_6m0j_max_escape.pdb
        results/pdb_outputs/LY-CoV555+LY-CoV016_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/LY-CoV555+LY-CoV016_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for LY-CoV555_400 to:
        results/pdb_outputs/LY-CoV555_400_6m0j_total_escape.pdb
        results/pdb_outputs/LY-CoV555_400_6m0j_max_escape.pdb
        results/pdb_outputs/LY-CoV555_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/LY-CoV555_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for REGN10933+REGN10987_400 to:
        results/pdb_outputs/REGN10933+REGN10987_400_6m0j_total_escape.pdb
        results/pdb_outputs/REGN10933+REGN10987_400_6m0j_max_escape.pdb
        results/pdb_outputs/REGN10933+REGN10987_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/REGN10933+REGN10987_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for REGN10933_400 to:
        results/pdb_outputs/REGN10933_400_6m0j_total_escape.pdb
        results/pdb_outputs/REGN10933_400_6m0j_max_escape.pdb
        results/pdb_outputs/REGN10933_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/REGN10933_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for REGN10987_400 to:
        results/pdb_outputs/REGN10987_400_6m0j_total_escape.pdb
        results/pdb_outputs/REGN10987_400_6m0j_max_escape.pdb
        results/pdb_outputs/REGN10987_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/REGN10987_400_6m0j_norm_max_escape.pdb
    
    Making PDB mappings for 6xcm to data/pdbs/6xcm.pdb
    Making mappings for 1 conditions.
    Mapping to the following chains: A, B, C
      Writing B-factor re-assigned PDBs for C105_400 to:
        results/pdb_outputs/C105_400_6xcm_total_escape.pdb
        results/pdb_outputs/C105_400_6xcm_max_escape.pdb
        results/pdb_outputs/C105_400_6xcm_norm_total_escape.pdb
        results/pdb_outputs/C105_400_6xcm_norm_max_escape.pdb
    
    Making PDB mappings for 6xcn to data/pdbs/6xcn.pdb
    Making mappings for 1 conditions.
    Mapping to the following chains: A, C, E
      Writing B-factor re-assigned PDBs for C105_400 to:
        results/pdb_outputs/C105_400_6xcn_total_escape.pdb
        results/pdb_outputs/C105_400_6xcn_max_escape.pdb
        results/pdb_outputs/C105_400_6xcn_norm_total_escape.pdb
        results/pdb_outputs/C105_400_6xcn_norm_max_escape.pdb
    
    Making PDB mappings for 7k90 to data/pdbs/7k90.pdb
    Making mappings for 1 conditions.
    Mapping to the following chains: A, B, C
      Writing B-factor re-assigned PDBs for C144_400 to:
        results/pdb_outputs/C144_400_7k90_total_escape.pdb
        results/pdb_outputs/C144_400_7k90_max_escape.pdb
        results/pdb_outputs/C144_400_7k90_norm_total_escape.pdb
        results/pdb_outputs/C144_400_7k90_norm_max_escape.pdb
    
    Making PDB mappings for 7k8s to data/pdbs/7k8s.pdb
    Making mappings for 1 conditions.
    Mapping to the following chains: A, B, C
      Writing B-factor re-assigned PDBs for C002_400 to:
        results/pdb_outputs/C002_400_7k8s_total_escape.pdb
        results/pdb_outputs/C002_400_7k8s_max_escape.pdb
        results/pdb_outputs/C002_400_7k8s_norm_total_escape.pdb
        results/pdb_outputs/C002_400_7k8s_norm_max_escape.pdb
    
    Making PDB mappings for 7k8t to data/pdbs/7k8t.pdb
    Making mappings for 1 conditions.
    Mapping to the following chains: A, B, C
      Writing B-factor re-assigned PDBs for C002_400 to:
        results/pdb_outputs/C002_400_7k8t_total_escape.pdb
        results/pdb_outputs/C002_400_7k8t_max_escape.pdb
        results/pdb_outputs/C002_400_7k8t_norm_total_escape.pdb
        results/pdb_outputs/C002_400_7k8t_norm_max_escape.pdb
    
    Making PDB mappings for 7k8x to data/pdbs/7k8x.pdb
    Making mappings for 1 conditions.
    Mapping to the following chains: A, B, C
      Writing B-factor re-assigned PDBs for C121_400 to:
        results/pdb_outputs/C121_400_7k8x_total_escape.pdb
        results/pdb_outputs/C121_400_7k8x_max_escape.pdb
        results/pdb_outputs/C121_400_7k8x_norm_total_escape.pdb
        results/pdb_outputs/C121_400_7k8x_norm_max_escape.pdb
    
    Making PDB mappings for 7k8y to data/pdbs/7k8y.pdb
    Making mappings for 1 conditions.
    Mapping to the following chains: B, D, E
      Writing B-factor re-assigned PDBs for C121_400 to:
        results/pdb_outputs/C121_400_7k8y_total_escape.pdb
        results/pdb_outputs/C121_400_7k8y_max_escape.pdb
        results/pdb_outputs/C121_400_7k8y_norm_total_escape.pdb
        results/pdb_outputs/C121_400_7k8y_norm_max_escape.pdb
    
    Making PDB mappings for 7k8z to data/pdbs/7k8z.pdb
    Making mappings for 1 conditions.
    Mapping to the following chains: A, B, C
      Writing B-factor re-assigned PDBs for C135_400 to:
        results/pdb_outputs/C135_400_7k8z_total_escape.pdb
        results/pdb_outputs/C135_400_7k8z_max_escape.pdb
        results/pdb_outputs/C135_400_7k8z_norm_total_escape.pdb
        results/pdb_outputs/C135_400_7k8z_norm_max_escape.pdb
    
    Making PDB mappings for 7k8v to data/pdbs/7k8v.pdb
    Making mappings for 1 conditions.
    Mapping to the following chains: A, B, C
      Writing B-factor re-assigned PDBs for C110_400 to:
        results/pdb_outputs/C110_400_7k8v_total_escape.pdb
        results/pdb_outputs/C110_400_7k8v_max_escape.pdb
        results/pdb_outputs/C110_400_7k8v_norm_total_escape.pdb
        results/pdb_outputs/C110_400_7k8v_norm_max_escape.pdb
    
    Making PDB mappings for 7c01 to data/pdbs/7c01_single.pdb
    Making mappings for 1 conditions.
    Mapping to the following chains: A
      Writing B-factor re-assigned PDBs for CB6_400 to:
        results/pdb_outputs/CB6_400_7c01_total_escape.pdb
        results/pdb_outputs/CB6_400_7c01_max_escape.pdb
        results/pdb_outputs/CB6_400_7c01_norm_total_escape.pdb
        results/pdb_outputs/CB6_400_7c01_norm_max_escape.pdb

