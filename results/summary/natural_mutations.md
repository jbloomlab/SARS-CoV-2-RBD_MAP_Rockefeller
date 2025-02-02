# Analyze naturally occurring mutations at sites of strong escape
This Python Jupyter notebook sees how many naturally occuring mutations are observed at each site of strong escape

## Set up analysis
Import Python modules:


```python
import collections
import copy
import math
import os

import dms_variants.utils

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import pandas as pd

from plotnine import *

import yaml
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Read escape profiles config, which tells which sets to make plots for:


```python
with open(config['escape_profiles_config']) as f:
    escape_profiles_config = yaml.safe_load(f)
```

Create output directory:


```python
os.makedirs(config['gisaid_mutations_dir'], exist_ok=True)
```

Read counts of naturally ocurring mutations:


```python
print(f"Reading mutation counts from {config['gisaid_mutation_counts']}")

mut_counts = pd.read_csv(config['gisaid_mutation_counts'])
```

    Reading mutation counts from results/GISAID_mutations/mutation_counts.csv


Read sites of "strong escape" from all antibodies / sera:


```python
print(f"Reading sites of strong escape from {config['strong_escape_sites']}")

strong_sites = pd.read_csv(config['strong_escape_sites'])
```

    Reading sites of strong escape from results/escape_profiles/strong_escape_sites.csv


Read escape fractions for all antibodies / sera:


```python
print(f"Reading escape fractions from {config['escape_fracs']}")

escape_fracs = (
    pd.read_csv(config['escape_fracs'])
    .query('library == "average"')
    .drop(columns='site')
    .rename(columns={'mutation': 'mutant',
                     'label_site': 'site'})
    [['condition', 'site', 'wildtype', 'mutant', config['mut_metric'], config['site_metric']]]
    )

escape_fracs
```

    Reading escape fractions from results/escape_scores/escape_fracs.csv





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>condition</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutant</th>
      <th>mut_escape_frac_epistasis_model</th>
      <th>site_total_escape_frac_epistasis_model</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>0.002020</td>
      <td>0.04926</td>
    </tr>
    <tr>
      <th>1</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>0.005616</td>
      <td>0.04926</td>
    </tr>
    <tr>
      <th>2</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>0.002535</td>
      <td>0.04926</td>
    </tr>
    <tr>
      <th>3</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>F</td>
      <td>0.003032</td>
      <td>0.04926</td>
    </tr>
    <tr>
      <th>4</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>G</td>
      <td>0.003113</td>
      <td>0.04926</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>103178</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>R</td>
      <td>0.002448</td>
      <td>0.04056</td>
    </tr>
    <tr>
      <th>103179</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>S</td>
      <td>0.002555</td>
      <td>0.04056</td>
    </tr>
    <tr>
      <th>103180</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>V</td>
      <td>0.002227</td>
      <td>0.04056</td>
    </tr>
    <tr>
      <th>103181</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>W</td>
      <td>0.001916</td>
      <td>0.04056</td>
    </tr>
    <tr>
      <th>103182</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>Y</td>
      <td>0.002564</td>
      <td>0.04056</td>
    </tr>
  </tbody>
</table>
<p>103183 rows × 6 columns</p>
</div>



## Counts of mutations at sites of escape
Get counts of naturally occurring mutations at sites of escape, along with the actual escape values:

First get mutation-level counts:


```python
mutcounts_strong_sites = (
    strong_sites[['condition', 'threshold', 'site']]
    .merge(mut_counts, how='inner', on='site')
    .merge(escape_fracs[['condition', 'site', 'wildtype', config['site_metric']]].drop_duplicates(),
           on=['condition', 'site', 'wildtype'],
           validate='many_to_one')
    .assign(mutation=lambda x: x['wildtype'] + x['site'].astype(str) + x['mutant'])
    .sort_values('count', ascending=False)
    )
```

Now get site-level counts (aggregating all mutations at a site):


```python
sitecounts_strong_sites = (
    mutcounts_strong_sites
    .assign(mut_count=lambda x: x['mutation'] + ' (' + x['count'].astype(str) + ')')
    .groupby(['condition', 'threshold', 'site', 'wildtype', config['site_metric']])
    .aggregate({'count': 'sum', 'mut_count': ', '.join})
    .rename(columns={'mut_count': 'counts_by_mutation'})
    .reset_index()
    .sort_values('count', ascending=False)
    )

print(f"Here are first few lines showing the most frequently mutated sites of escape:")
display(HTML(sitecounts_strong_sites.head(n=20).to_html(index=False)))
```

    Here are first few lines showing the most frequently mutated sites of escape:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>condition</th>
      <th>threshold</th>
      <th>site</th>
      <th>wildtype</th>
      <th>site_total_escape_frac_epistasis_model</th>
      <th>count</th>
      <th>counts_by_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>COV2-2499_400</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>0.4765</td>
      <td>563243</td>
      <td>N501Y (560239), N501T (2938), N501S (31), N501I (19), N501H (13), N501C (1), N501F (1), N501K (1)</td>
    </tr>
    <tr>
      <td>25_d18_500</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>0.6856</td>
      <td>563243</td>
      <td>N501Y (560239), N501T (2938), N501S (31), N501I (19), N501H (13), N501K (1), N501F (1), N501C (1)</td>
    </tr>
    <tr>
      <td>22C_d28_200</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>1.2370</td>
      <td>563243</td>
      <td>N501Y (560239), N501T (2938), N501S (31), N501I (19), N501H (13), N501K (1), N501F (1), N501C (1)</td>
    </tr>
    <tr>
      <td>COV-021_500</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>0.3541</td>
      <td>563243</td>
      <td>N501Y (560239), N501T (2938), N501S (31), N501I (19), N501H (13), N501K (1), N501F (1), N501C (1)</td>
    </tr>
    <tr>
      <td>COV-021_500</td>
      <td>sensitive_max_mut</td>
      <td>501</td>
      <td>N</td>
      <td>0.3541</td>
      <td>563243</td>
      <td>N501Y (560239), N501T (2938), N501S (31), N501I (19), N501H (13), N501K (1), N501F (1), N501C (1)</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>sensitive_max_mut</td>
      <td>501</td>
      <td>N</td>
      <td>1.5790</td>
      <td>563243</td>
      <td>N501Y (560239), N501T (2938), N501S (31), N501I (19), N501H (13), N501K (1), N501F (1), N501C (1)</td>
    </tr>
    <tr>
      <td>25C_d115_80</td>
      <td>sensitive_max_mut</td>
      <td>501</td>
      <td>N</td>
      <td>0.9038</td>
      <td>563243</td>
      <td>N501Y (560239), N501T (2938), N501S (31), N501I (19), N501H (13), N501K (1), N501F (1), N501C (1)</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>1.5790</td>
      <td>563243</td>
      <td>N501Y (560239), N501T (2938), N501S (31), N501I (19), N501H (13), N501K (1), N501F (1), N501C (1)</td>
    </tr>
    <tr>
      <td>25C_d115_80</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>0.9038</td>
      <td>563243</td>
      <td>N501Y (560239), N501T (2938), N501S (31), N501I (19), N501H (13), N501F (1), N501C (1), N501K (1)</td>
    </tr>
    <tr>
      <td>22C_d28_200</td>
      <td>sensitive_max_mut</td>
      <td>501</td>
      <td>N</td>
      <td>1.2370</td>
      <td>563243</td>
      <td>N501Y (560239), N501T (2938), N501S (31), N501I (19), N501H (13), N501K (1), N501F (1), N501C (1)</td>
    </tr>
    <tr>
      <td>COV2-2499_400</td>
      <td>sensitive_max_mut</td>
      <td>501</td>
      <td>N</td>
      <td>0.4765</td>
      <td>563243</td>
      <td>N501Y (560239), N501T (2938), N501S (31), N501I (19), N501H (13), N501K (1), N501F (1), N501C (1)</td>
    </tr>
    <tr>
      <td>25_d18_500</td>
      <td>sensitive_max_mut</td>
      <td>501</td>
      <td>N</td>
      <td>0.6856</td>
      <td>563243</td>
      <td>N501Y (560239), N501T (2938), N501S (31), N501I (19), N501H (13), N501K (1), N501F (1), N501C (1)</td>
    </tr>
    <tr>
      <td>C144_400</td>
      <td>sensitive</td>
      <td>484</td>
      <td>E</td>
      <td>18.4400</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>12C_d61_160</td>
      <td>sensitive_max_mut</td>
      <td>484</td>
      <td>E</td>
      <td>1.5620</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>25C_d115_80</td>
      <td>sensitive</td>
      <td>484</td>
      <td>E</td>
      <td>4.9770</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>13_d15_200</td>
      <td>sensitive_max_mut</td>
      <td>484</td>
      <td>E</td>
      <td>0.5454</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>25_d18_500</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>2.7910</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>25C_d48_200</td>
      <td>sensitive_max_mut</td>
      <td>484</td>
      <td>E</td>
      <td>4.1490</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>COV-021_500</td>
      <td>sensitive</td>
      <td>484</td>
      <td>E</td>
      <td>2.3400</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>COV-021_500</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>2.3400</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
  </tbody>
</table>


## Plot sites of escape with natural variation
We perform analyses on all subsets in the escape profiles config for which this is specified:


```python
for name, specs in escape_profiles_config.items():
    if 'analyze_natural_mutations' not in specs or not specs['analyze_natural_mutations']:
        continue
    print(f"\nAnalyzing natural mutations for {name}")
    
    conditions = specs['conditions']
    
    threshold = specs['plot_auto_identified_sites']
    if threshold not in sitecounts_strong_sites['threshold'].unique():
        raise ValueError(f"invalid threshold {threshold} for {name}")
    
    # get count for conditions of interest for this subset
    df = (sitecounts_strong_sites
          .query('condition in @conditions')
          .query('threshold == @threshold')
          .assign(condition=lambda x: x['condition'].map(conditions))
          .drop(columns=config['site_metric'])
          )
    countsfile = os.path.join(config['gisaid_mutations_dir'], f"{name}_mutation_counts.csv")
    print(f"Writing counts of mutations at sites of strong escape to {countsfile}. First few lines:")
    display(HTML(df.head(n=10).to_html(index=False)))
    df.to_csv(countsfile, index=False)
    
    # make plot showing escape sites with more than mincounts mutations
    if 'natural_mutations_mincounts' in specs:
        mincounts = specs['natural_mutations_mincounts']
    else:
        mincounts = 5
    plotsfile = os.path.join(config['gisaid_mutations_dir'], f"{name}_mutation_counts.pdf")
    print('Plotting which antibodies / sera are escaped by mutations at all sites of '
          f"escape with at least {mincounts} mutation counts and saving to {plotsfile}.")
    plot_df = (
        # data frame with all combinations of conditions and sites
        pd.DataFrame.from_records([(condition, site) for condition in conditions.values()
                                   for site in df['site'].unique()],
                                  columns=['condition', 'site'])
        # annotate sites of escape
        .merge(df.assign(escape=lambda x: x['count'] >= mincounts)[['condition', 'site', 'escape']],
               how='left',
               validate='one_to_one',
               on=['condition', 'site'])
        .assign(escape=lambda x: x['escape'].fillna(False))
        # add wildtype and counts of mutations at each site
        .merge(sitecounts_strong_sites[['site', 'wildtype', 'count']].drop_duplicates(),
               how='left',
               validate='many_to_one',
               on='site')
        # get only sites with sufficient mutation counts
        .query('count > @mincounts')
        # only get sites where at least one antibody escapes
        .assign(n_escape=lambda x: x.groupby('site')['escape'].transform('sum'))
        .query('n_escape > 0')
        # order conditions, and order sites by count after making nice label
        .assign(site_label=lambda x: x['wildtype'] + x['site'].astype(str) + ' (' + x['count'].astype(str) + ')')
        .sort_values('count')
        .assign(condition=lambda x: pd.Categorical(x['condition'], list(conditions.values()), ordered=True),
                site_label=lambda x: pd.Categorical(x['site_label'], x['site_label'].unique(), ordered=True)
                )
        )
    p = (ggplot(plot_df) +
         aes('condition', 'site_label', fill='escape') +
         geom_tile(color='black', size=0.3) +
         theme(axis_text_x=element_text(angle=90),
               figure_size=(0.3 * plot_df['condition'].nunique(), 0.3 * plot_df['site_label'].nunique()),
               panel_background=element_blank(),
               ) +
         xlab('') +
         ylab('') +
         scale_fill_manual(values=['white', 'dimgray'])
         )
    p.save(plotsfile, verbose=False)
    fig = p.draw()
    display(fig)
    plt.close(fig)
```

    
    Analyzing natural mutations for Nussenzweig_serum
    Writing counts of mutations at sites of strong escape to results/GISAID_mutations/Nussenzweig_serum_mutation_counts.csv. First few lines:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>condition</th>
      <th>threshold</th>
      <th>site</th>
      <th>wildtype</th>
      <th>count</th>
      <th>counts_by_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>COV-021</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>COV-072</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>COV-057</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>COV-047</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>COV-047</td>
      <td>default</td>
      <td>490</td>
      <td>F</td>
      <td>2195</td>
      <td>F490S (1892), F490L (217), F490V (38), F490Y (38), F490I (8), F490R (1), F490P (1)</td>
    </tr>
    <tr>
      <td>COV-057</td>
      <td>default</td>
      <td>490</td>
      <td>F</td>
      <td>2195</td>
      <td>F490S (1892), F490L (217), F490V (38), F490Y (38), F490I (8), F490R (1), F490P (1)</td>
    </tr>
    <tr>
      <td>COV-057</td>
      <td>default</td>
      <td>446</td>
      <td>G</td>
      <td>471</td>
      <td>G446V (382), G446S (73), G446D (8), G446A (7), G446R (1)</td>
    </tr>
    <tr>
      <td>COV-047</td>
      <td>default</td>
      <td>455</td>
      <td>L</td>
      <td>405</td>
      <td>L455F (394), L455S (9), L455V (2)</td>
    </tr>
    <tr>
      <td>COV-021</td>
      <td>default</td>
      <td>455</td>
      <td>L</td>
      <td>405</td>
      <td>L455F (394), L455S (9), L455V (2)</td>
    </tr>
    <tr>
      <td>COV-057</td>
      <td>default</td>
      <td>450</td>
      <td>N</td>
      <td>318</td>
      <td>N450K (295), N450D (18), N450S (3), N450I (1), N450H (1)</td>
    </tr>
  </tbody>
</table>


    Plotting which antibodies / sera are escaped by mutations at all sites of escape with at least 5 mutation counts and saving to results/GISAID_mutations/Nussenzweig_serum_mutation_counts.pdf.



    
![png](natural_mutations_files/natural_mutations_20_3.png)
    


    
    Analyzing natural mutations for all_class1_abs
    Writing counts of mutations at sites of strong escape to results/GISAID_mutations/all_class1_abs_mutation_counts.csv. First few lines:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>condition</th>
      <th>threshold</th>
      <th>site</th>
      <th>wildtype</th>
      <th>count</th>
      <th>counts_by_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>COV2-2832</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>REGN10933</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>C105</td>
      <td>default</td>
      <td>417</td>
      <td>K</td>
      <td>22447</td>
      <td>K417N (12354), K417T (10070), K417R (18), K417E (3), K417S (1), K417M (1)</td>
    </tr>
    <tr>
      <td>LY-CoV016</td>
      <td>default</td>
      <td>417</td>
      <td>K</td>
      <td>22447</td>
      <td>K417N (12354), K417T (10070), K417R (18), K417E (3), K417M (1), K417S (1)</td>
    </tr>
    <tr>
      <td>REGN10933</td>
      <td>default</td>
      <td>417</td>
      <td>K</td>
      <td>22447</td>
      <td>K417N (12354), K417T (10070), K417R (18), K417E (3), K417S (1), K417M (1)</td>
    </tr>
    <tr>
      <td>REGN10933</td>
      <td>default</td>
      <td>453</td>
      <td>Y</td>
      <td>1206</td>
      <td>Y453F (1200), Y453H (4), Y453C (1), Y453S (1)</td>
    </tr>
    <tr>
      <td>C105</td>
      <td>default</td>
      <td>453</td>
      <td>Y</td>
      <td>1206</td>
      <td>Y453F (1200), Y453H (4), Y453S (1), Y453C (1)</td>
    </tr>
    <tr>
      <td>C105</td>
      <td>default</td>
      <td>475</td>
      <td>A</td>
      <td>425</td>
      <td>A475V (300), A475S (98), A475T (24), A475D (2), A475P (1)</td>
    </tr>
    <tr>
      <td>REGN10933</td>
      <td>default</td>
      <td>475</td>
      <td>A</td>
      <td>425</td>
      <td>A475V (300), A475S (98), A475T (24), A475D (2), A475P (1)</td>
    </tr>
    <tr>
      <td>COV2-2165</td>
      <td>default</td>
      <td>475</td>
      <td>A</td>
      <td>425</td>
      <td>A475V (300), A475S (98), A475T (24), A475D (2), A475P (1)</td>
    </tr>
  </tbody>
</table>


    Plotting which antibodies / sera are escaped by mutations at all sites of escape with at least 5 mutation counts and saving to results/GISAID_mutations/all_class1_abs_mutation_counts.pdf.



    
![png](natural_mutations_files/natural_mutations_20_7.png)
    


    
    Analyzing natural mutations for all_class2_abs
    Writing counts of mutations at sites of strong escape to results/GISAID_mutations/all_class2_abs_mutation_counts.csv. First few lines:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>condition</th>
      <th>threshold</th>
      <th>site</th>
      <th>wildtype</th>
      <th>count</th>
      <th>counts_by_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>C144</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>COV2-2479</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>COV2-2050</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>C002</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>C121</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>LY-CoV555</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>LY-CoV555</td>
      <td>default</td>
      <td>452</td>
      <td>L</td>
      <td>50704</td>
      <td>L452R (49525), L452Q (833), L452M (340), L452V (3), L452P (2), L452F (1)</td>
    </tr>
    <tr>
      <td>C002</td>
      <td>default</td>
      <td>452</td>
      <td>L</td>
      <td>50704</td>
      <td>L452R (49525), L452Q (833), L452M (340), L452V (3), L452P (2), L452F (1)</td>
    </tr>
    <tr>
      <td>LY-CoV555</td>
      <td>default</td>
      <td>494</td>
      <td>S</td>
      <td>6550</td>
      <td>S494P (6469), S494L (60), S494A (14), S494T (4), S494R (3)</td>
    </tr>
    <tr>
      <td>C121</td>
      <td>default</td>
      <td>490</td>
      <td>F</td>
      <td>2195</td>
      <td>F490S (1892), F490L (217), F490V (38), F490Y (38), F490I (8), F490P (1), F490R (1)</td>
    </tr>
  </tbody>
</table>


    Plotting which antibodies / sera are escaped by mutations at all sites of escape with at least 5 mutation counts and saving to results/GISAID_mutations/all_class2_abs_mutation_counts.pdf.



    
![png](natural_mutations_files/natural_mutations_20_11.png)
    


    
    Analyzing natural mutations for all_class3_abs
    Writing counts of mutations at sites of strong escape to results/GISAID_mutations/all_class3_abs_mutation_counts.csv. First few lines:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>condition</th>
      <th>threshold</th>
      <th>site</th>
      <th>wildtype</th>
      <th>count</th>
      <th>counts_by_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>C110</td>
      <td>default</td>
      <td>452</td>
      <td>L</td>
      <td>50704</td>
      <td>L452R (49525), L452Q (833), L452M (340), L452V (3), L452P (2), L452F (1)</td>
    </tr>
    <tr>
      <td>REGN10987</td>
      <td>default</td>
      <td>439</td>
      <td>N</td>
      <td>19872</td>
      <td>N439K (19861), N439I (6), N439S (2), N439Y (1), N439T (1), N439D (1)</td>
    </tr>
    <tr>
      <td>COV2-2130</td>
      <td>default</td>
      <td>494</td>
      <td>S</td>
      <td>6550</td>
      <td>S494P (6469), S494L (60), S494A (14), S494T (4), S494R (3)</td>
    </tr>
    <tr>
      <td>C110</td>
      <td>default</td>
      <td>494</td>
      <td>S</td>
      <td>6550</td>
      <td>S494P (6469), S494L (60), S494A (14), S494T (4), S494R (3)</td>
    </tr>
    <tr>
      <td>C110</td>
      <td>default</td>
      <td>490</td>
      <td>F</td>
      <td>2195</td>
      <td>F490S (1892), F490L (217), F490V (38), F490Y (38), F490I (8), F490P (1), F490R (1)</td>
    </tr>
    <tr>
      <td>C135</td>
      <td>default</td>
      <td>440</td>
      <td>N</td>
      <td>1877</td>
      <td>N440K (1783), N440Y (37), N440T (21), N440S (14), N440D (11), N440I (7), N440F (2), N440H (1), N440R (1)</td>
    </tr>
    <tr>
      <td>REGN10987</td>
      <td>default</td>
      <td>440</td>
      <td>N</td>
      <td>1877</td>
      <td>N440K (1783), N440Y (37), N440T (21), N440S (14), N440D (11), N440I (7), N440F (2), N440R (1), N440H (1)</td>
    </tr>
    <tr>
      <td>COV2-2130</td>
      <td>default</td>
      <td>346</td>
      <td>R</td>
      <td>1019</td>
      <td>R346K (573), R346S (331), R346I (63), R346G (35), R346T (17)</td>
    </tr>
    <tr>
      <td>C110</td>
      <td>default</td>
      <td>346</td>
      <td>R</td>
      <td>1019</td>
      <td>R346K (573), R346S (331), R346I (63), R346G (35), R346T (17)</td>
    </tr>
    <tr>
      <td>C135</td>
      <td>default</td>
      <td>346</td>
      <td>R</td>
      <td>1019</td>
      <td>R346K (573), R346S (331), R346I (63), R346G (35), R346T (17)</td>
    </tr>
  </tbody>
</table>


    Plotting which antibodies / sera are escaped by mutations at all sites of escape with at least 5 mutation counts and saving to results/GISAID_mutations/all_class3_abs_mutation_counts.pdf.



    
![png](natural_mutations_files/natural_mutations_20_15.png)
    


    
    Analyzing natural mutations for all_class4_abs
    Writing counts of mutations at sites of strong escape to results/GISAID_mutations/all_class4_abs_mutation_counts.csv. First few lines:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>condition</th>
      <th>threshold</th>
      <th>site</th>
      <th>wildtype</th>
      <th>count</th>
      <th>counts_by_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>COV2-2094</td>
      <td>default</td>
      <td>417</td>
      <td>K</td>
      <td>22447</td>
      <td>K417N (12354), K417T (10070), K417R (18), K417E (3), K417S (1), K417M (1)</td>
    </tr>
    <tr>
      <td>COV2-2082</td>
      <td>default</td>
      <td>417</td>
      <td>K</td>
      <td>22447</td>
      <td>K417N (12354), K417T (10070), K417R (18), K417E (3), K417S (1), K417M (1)</td>
    </tr>
    <tr>
      <td>COV2-2677</td>
      <td>default</td>
      <td>384</td>
      <td>P</td>
      <td>1927</td>
      <td>P384L (1325), P384S (560), P384R (31), P384H (6), P384T (3), P384A (2)</td>
    </tr>
    <tr>
      <td>COV2-2094</td>
      <td>default</td>
      <td>408</td>
      <td>R</td>
      <td>435</td>
      <td>R408K (218), R408I (138), R408G (63), R408T (13), R408S (3)</td>
    </tr>
    <tr>
      <td>COV2-2082</td>
      <td>default</td>
      <td>408</td>
      <td>R</td>
      <td>435</td>
      <td>R408K (218), R408I (138), R408G (63), R408T (13), R408S (3)</td>
    </tr>
    <tr>
      <td>COV2-2677</td>
      <td>default</td>
      <td>370</td>
      <td>N</td>
      <td>259</td>
      <td>N370S (123), N370K (83), N370H (46), N370D (7)</td>
    </tr>
    <tr>
      <td>COV2-2082</td>
      <td>default</td>
      <td>376</td>
      <td>T</td>
      <td>172</td>
      <td>T376I (150), T376S (13), T376N (9)</td>
    </tr>
    <tr>
      <td>COV2-2094</td>
      <td>default</td>
      <td>376</td>
      <td>T</td>
      <td>172</td>
      <td>T376I (150), T376S (13), T376N (9)</td>
    </tr>
    <tr>
      <td>COV2-2677</td>
      <td>default</td>
      <td>372</td>
      <td>A</td>
      <td>108</td>
      <td>A372V (82), A372T (10), A372P (7), A372S (6), A372L (2), A372G (1)</td>
    </tr>
    <tr>
      <td>COV2-2677</td>
      <td>default</td>
      <td>378</td>
      <td>K</td>
      <td>99</td>
      <td>K378N (81), K378R (12), K378M (5), K378E (1)</td>
    </tr>
  </tbody>
</table>


    Plotting which antibodies / sera are escaped by mutations at all sites of escape with at least 5 mutation counts and saving to results/GISAID_mutations/all_class4_abs_mutation_counts.pdf.



    
![png](natural_mutations_files/natural_mutations_20_19.png)
    


    
    Analyzing natural mutations for HAARVI_sera
    Writing counts of mutations at sites of strong escape to results/GISAID_mutations/HAARVI_sera_mutation_counts.csv. First few lines:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>condition</th>
      <th>threshold</th>
      <th>site</th>
      <th>wildtype</th>
      <th>count</th>
      <th>counts_by_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>subject G (day 18)</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>subject F (day 48)</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>subject I (day 102)</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>subject B (day 113)</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>subject E (day 104)</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>subject A (day 21)</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>subject B (day 26)</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>subject E (day 28)</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>subject I (day 26)</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
    <tr>
      <td>subject C (day 104)</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>52154</td>
      <td>E484K (50055), E484Q (1978), E484G (51), E484A (36), E484D (28), E484R (4), E484V (2)</td>
    </tr>
  </tbody>
</table>


    Plotting which antibodies / sera are escaped by mutations at all sites of escape with at least 5 mutation counts and saving to results/GISAID_mutations/HAARVI_sera_mutation_counts.pdf.



    
![png](natural_mutations_files/natural_mutations_20_23.png)
    


## Plot correlation between escape and natural frequency
First aggregate frequency of mutations and escape fractions:


```python
escape_and_freq = (
    escape_fracs
    .rename(columns={config['mut_metric']: 'mut_escape',
                     config['site_metric']: 'tot_site_escape'})
    .assign(max_site_escape=lambda x: x.groupby(['condition', 'site'])['mut_escape'].transform('max'),
            mean_site_escape=lambda x: x.groupby(['condition', 'site'])['mut_escape'].transform('mean'))
    .merge(mut_counts[['site', 'wildtype', 'mutant', 'frequency']]
                     .rename(columns={'frequency': 'mut_freq'}),
           on=['site', 'wildtype', 'mutant'],
           how='left', validate='many_to_one')
    .assign(mut_freq=lambda x: x['mut_freq'].fillna(0),
            site_freq=lambda x: x.groupby(['condition', 'site'])['mut_freq'].transform('sum'),
            mutation=lambda x: x['wildtype'] + x['site'].astype(str) + x['mutant'],
            )
    )

display(HTML(escape_and_freq.head().to_html()))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>condition</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutant</th>
      <th>mut_escape</th>
      <th>tot_site_escape</th>
      <th>max_site_escape</th>
      <th>mean_site_escape</th>
      <th>mut_freq</th>
      <th>site_freq</th>
      <th>mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>0.002020</td>
      <td>0.04926</td>
      <td>0.007632</td>
      <td>0.003079</td>
      <td>0.000000e+00</td>
      <td>0.000032</td>
      <td>N331A</td>
    </tr>
    <tr>
      <th>1</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>0.005616</td>
      <td>0.04926</td>
      <td>0.007632</td>
      <td>0.003079</td>
      <td>7.531181e-07</td>
      <td>0.000032</td>
      <td>N331D</td>
    </tr>
    <tr>
      <th>2</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>0.002535</td>
      <td>0.04926</td>
      <td>0.007632</td>
      <td>0.003079</td>
      <td>0.000000e+00</td>
      <td>0.000032</td>
      <td>N331E</td>
    </tr>
    <tr>
      <th>3</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>F</td>
      <td>0.003032</td>
      <td>0.04926</td>
      <td>0.007632</td>
      <td>0.003079</td>
      <td>0.000000e+00</td>
      <td>0.000032</td>
      <td>N331F</td>
    </tr>
    <tr>
      <th>4</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>G</td>
      <td>0.003113</td>
      <td>0.04926</td>
      <td>0.007632</td>
      <td>0.003079</td>
      <td>0.000000e+00</td>
      <td>0.000032</td>
      <td>N331G</td>
    </tr>
  </tbody>
</table>


Now make plots.
Note that you can configure below exactly what variables you want to plot (mutation frequency, mutation escape, site escape, etc):


```python
# all the parameters below have the indicated defaults, but can be set in `escape_profiles_config`
# via analyze_natural_mutations_specs
default_analysis_specs = {
    'maxcol': 5,  # maximum columns in plot
    'minfreq': 1e-5,  # collapse any natural frequencies < this
    'freq': 'site_freq',  # type of frequency to plot: mut_freq or site_freq
    'escape': 'mean_site_escape',  # type of escape to plot: mean_site_escape, mut_escape, max_site_escape, tot_site_escape
    'xlabel': 'frequency of mutations at site',
    'ylabel': 'mean escape at site',
    'label_minfreq': 5e-5,  # label points with frequency >= this and...
    'label_minescape': 0.1,  # label points with escape >= this
    'also_label': [],  # also label any points (sites or mutations) listed here
    'label_font_size': 6,  # font size for labeling points
    'default_color': 'black',  # color for points not otherwise specified
    'default_alpha': 0.5,  # default alpha if not specified
    'set_point_color': {},  # set color; key by site / mutation, value is color
    'set_point_alpha': {},  # set alpha; key by site / mutations, value is alpha
    'plot_average_only': False,  # only plot average of conditions, not individual ones
    }
label_minfreq = 5e-5  # label points with frequency >= this
label_minescape = 0.05  # label points with escape >= this

for name, specs in escape_profiles_config.items():
    if 'analyze_natural_mutations' not in specs or not specs['analyze_natural_mutations']:
        continue
    print(f"\nAnalyzing natural mutations for {name}")
    
    analysis_specs = copy.deepcopy(default_analysis_specs)
    if 'analyze_natural_mutations_specs' in specs:
        for key, val in specs['analyze_natural_mutations_specs'].items():
            analysis_specs[key] = val
    
    conditions = specs['conditions']
    
    if 'site' in analysis_specs['freq'] and 'site' in analysis_specs['escape']:
        ptlabel = 'site'
    else:
        ptlabel = 'mutation'
    
    df = (escape_and_freq
          .query('condition in @conditions')
          .assign(condition=lambda x: x['condition'].map(conditions))
          .assign(**{analysis_specs['freq']: lambda x: x[analysis_specs['freq']].clip(lower=analysis_specs['minfreq'])})
          [['condition', analysis_specs['escape'], analysis_specs['freq'], ptlabel]]
          .drop_duplicates()
          )

    assert len(conditions) == df['condition'].nunique()
    
    for avg_conditions in (False, True):
        
        if analysis_specs['plot_average_only'] and not avg_conditions:
            continue
        
        if avg_conditions:
            plot_df = df.groupby(ptlabel, as_index=False).aggregate({analysis_specs['freq']: 'mean',
                                                                     analysis_specs['escape']: 'mean'})
            nrow = ncol = 1
            plotfile = os.path.join(config['gisaid_mutations_dir'],
                                    f"{name}_escape_vs_freq_average.pdf")
            print(f"Plotting average across conditions and saving to {plotfile}")
        else:
            nrow = math.ceil(len(conditions) / analysis_specs['maxcol'])
            ncol = min(len(conditions), analysis_specs['maxcol'])
            plot_df = df.copy()
            # make condition categorical to maintain order 
            plot_df=plot_df.assign(condition=lambda x: pd.Categorical(x['condition'], 
                                                                      list(conditions.values()), 
                                                                      ordered=True)
                                  )
            plotfile = os.path.join(config['gisaid_mutations_dir'],
                                    f"{name}_escape_vs_freq_by-condition.pdf")
            print(f"Plotting each condition and saving to {plotfile}")
         
        # color points and set alpha
        set_point_color = collections.defaultdict(lambda: analysis_specs['default_color'])
        set_point_alpha = collections.defaultdict(lambda: analysis_specs['default_alpha'])
        for point, color in analysis_specs['set_point_color'].items():
            set_point_color[point] = color
        for point, alpha in analysis_specs['set_point_alpha'].items():
            set_point_alpha[point] = alpha
        plot_df['color'] = plot_df[ptlabel].map(set_point_color)
        plot_df['alpha'] = plot_df[ptlabel].map(set_point_alpha)
        # need to make color categorical to assign as aesthetic
        colors = plot_df['color'].unique()
        plot_df['color'] = pd.Categorical(plot_df['color'], colors, ordered=True)
            
        label_df = (plot_df
                    .assign(label=lambda x: x[ptlabel].isin(analysis_specs['also_label']))
                    .query(f"label or ({analysis_specs['freq']} >= {analysis_specs['label_minfreq']})")
                    .query(f"label or ({analysis_specs['escape']} >= {analysis_specs['label_minescape']})")
                    )
        
        maxfreq = plot_df[analysis_specs['freq']].max()
        assert analysis_specs['minfreq'] == 10**(int(math.log10(analysis_specs['minfreq'])))
        logxbreaks = list(range(int(math.log10(analysis_specs['minfreq'])), round(math.log10(maxfreq)) + 1, 1))
        xbreaks = [10**logx for logx in logxbreaks]
        xlabels = [f"$10^{{{logx}}}$" for logx in logxbreaks]
        xlabels[0] = f"$<{xlabels[0][1:]}"
        
        p = (ggplot(plot_df) +
             aes(analysis_specs['freq'], analysis_specs['escape'], color='color', alpha='alpha') +
             geom_point() +
             geom_text(data=label_df,
                       mapping=aes(label=ptlabel),
                       size=analysis_specs['label_font_size'],
                       adjust_text={'expand_points': (1.05, 1.2),
                                    'expand_text': (1.05, 1.2)},
                       ) +
             theme_classic() +
             theme(figure_size=(2.5 * ncol, 2.5 * nrow),
                   panel_spacing=0.3,
                   legend_position='none',
                   ) +
             scale_x_log10(name=analysis_specs['xlabel'],
                           breaks=xbreaks,
                           labels=xlabels,
                           expand=(0.07, 0)) +
             ylab(analysis_specs['ylabel']) +
             scale_color_manual(values=colors) +
             scale_alpha_continuous(limits=(0, 1), range=(0, 1))
             )
        if not avg_conditions:
            p = p + facet_wrap('~ condition', ncol=ncol, scales='free_y')
        p.save(plotfile, verbose=False)
        fig = p.draw()
        display(fig)
       # plt.close(fig)
```

    
    Analyzing natural mutations for Nussenzweig_serum
    Plotting average across conditions and saving to results/GISAID_mutations/Nussenzweig_serum_escape_vs_freq_average.pdf



    
![png](natural_mutations_files/natural_mutations_24_1.png)
    


    
    Analyzing natural mutations for all_class1_abs
    Plotting each condition and saving to results/GISAID_mutations/all_class1_abs_escape_vs_freq_by-condition.pdf


    /fh/fast/bloom_j/software/miniconda3/envs/SARS-CoV-2-RBD_MAP/lib/python3.7/site-packages/plotnine/facets/facet.py:552: PlotnineWarning: If you need more space for the x-axis tick text use ... + theme(subplots_adjust={'wspace': 0.25}). Choose an appropriate value for 'wspace'.
    /fh/fast/bloom_j/software/miniconda3/envs/SARS-CoV-2-RBD_MAP/lib/python3.7/site-packages/plotnine/facets/facet.py:552: PlotnineWarning: If you need more space for the x-axis tick text use ... + theme(subplots_adjust={'wspace': 0.25}). Choose an appropriate value for 'wspace'.



    
![png](natural_mutations_files/natural_mutations_24_4.png)
    


    Plotting average across conditions and saving to results/GISAID_mutations/all_class1_abs_escape_vs_freq_average.pdf



    
![png](natural_mutations_files/natural_mutations_24_6.png)
    


    
    Analyzing natural mutations for all_class2_abs
    Plotting each condition and saving to results/GISAID_mutations/all_class2_abs_escape_vs_freq_by-condition.pdf


    /fh/fast/bloom_j/software/miniconda3/envs/SARS-CoV-2-RBD_MAP/lib/python3.7/site-packages/plotnine/facets/facet.py:552: PlotnineWarning: If you need more space for the x-axis tick text use ... + theme(subplots_adjust={'wspace': 0.25}). Choose an appropriate value for 'wspace'.
    /fh/fast/bloom_j/software/miniconda3/envs/SARS-CoV-2-RBD_MAP/lib/python3.7/site-packages/plotnine/facets/facet.py:552: PlotnineWarning: If you need more space for the x-axis tick text use ... + theme(subplots_adjust={'wspace': 0.25}). Choose an appropriate value for 'wspace'.



    
![png](natural_mutations_files/natural_mutations_24_9.png)
    


    Plotting average across conditions and saving to results/GISAID_mutations/all_class2_abs_escape_vs_freq_average.pdf



    
![png](natural_mutations_files/natural_mutations_24_11.png)
    


    
    Analyzing natural mutations for all_class3_abs
    Plotting each condition and saving to results/GISAID_mutations/all_class3_abs_escape_vs_freq_by-condition.pdf


    /fh/fast/bloom_j/software/miniconda3/envs/SARS-CoV-2-RBD_MAP/lib/python3.7/site-packages/plotnine/facets/facet.py:552: PlotnineWarning: If you need more space for the x-axis tick text use ... + theme(subplots_adjust={'wspace': 0.25}). Choose an appropriate value for 'wspace'.
    /fh/fast/bloom_j/software/miniconda3/envs/SARS-CoV-2-RBD_MAP/lib/python3.7/site-packages/plotnine/facets/facet.py:552: PlotnineWarning: If you need more space for the x-axis tick text use ... + theme(subplots_adjust={'wspace': 0.25}). Choose an appropriate value for 'wspace'.



    
![png](natural_mutations_files/natural_mutations_24_14.png)
    


    Plotting average across conditions and saving to results/GISAID_mutations/all_class3_abs_escape_vs_freq_average.pdf



    
![png](natural_mutations_files/natural_mutations_24_16.png)
    


    
    Analyzing natural mutations for all_class4_abs
    Plotting each condition and saving to results/GISAID_mutations/all_class4_abs_escape_vs_freq_by-condition.pdf


    /fh/fast/bloom_j/software/miniconda3/envs/SARS-CoV-2-RBD_MAP/lib/python3.7/site-packages/plotnine/facets/facet.py:552: PlotnineWarning: If you need more space for the x-axis tick text use ... + theme(subplots_adjust={'wspace': 0.25}). Choose an appropriate value for 'wspace'.
    /fh/fast/bloom_j/software/miniconda3/envs/SARS-CoV-2-RBD_MAP/lib/python3.7/site-packages/plotnine/facets/facet.py:552: PlotnineWarning: If you need more space for the x-axis tick text use ... + theme(subplots_adjust={'wspace': 0.25}). Choose an appropriate value for 'wspace'.



    
![png](natural_mutations_files/natural_mutations_24_19.png)
    


    Plotting average across conditions and saving to results/GISAID_mutations/all_class4_abs_escape_vs_freq_average.pdf



    
![png](natural_mutations_files/natural_mutations_24_21.png)
    


    
    Analyzing natural mutations for HAARVI_sera
    Plotting average across conditions and saving to results/GISAID_mutations/HAARVI_sera_escape_vs_freq_average.pdf



    
![png](natural_mutations_files/natural_mutations_24_23.png)
    



    
![png](natural_mutations_files/natural_mutations_24_24.png)
    



    
![png](natural_mutations_files/natural_mutations_24_25.png)
    



    
![png](natural_mutations_files/natural_mutations_24_26.png)
    



    
![png](natural_mutations_files/natural_mutations_24_27.png)
    



    
![png](natural_mutations_files/natural_mutations_24_28.png)
    



    
![png](natural_mutations_files/natural_mutations_24_29.png)
    



    
![png](natural_mutations_files/natural_mutations_24_30.png)
    



    
![png](natural_mutations_files/natural_mutations_24_31.png)
    



    
![png](natural_mutations_files/natural_mutations_24_32.png)
    



    
![png](natural_mutations_files/natural_mutations_24_33.png)
    



```python

```
