# Align spike sequences in GISAID and count RBD mutations
This Python Jupyter notebook reads in a file of all spike sequences from GISAID, parses for "high-quality" sequences, builds a RBD alignment, and then makes a file that gives the counts of each mutation at each site.

## Set up analysis
Import Python modules:


```python
import io
import lzma
import os
import re
import subprocess

from Bio.Data.IUPACData import protein_letters
import Bio.SeqIO

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

Create output directory:


```python
os.makedirs(config['gisaid_mutations_dir'], exist_ok=True)
```

## Parse full-length human human spike sequences

Read the spikes from the file downloaded from GISAID:


```python
print(f"Reading GISAID spikes in {config['gisaid_spikes']}")
# file is `xz` compressed
with lzma.open(config['gisaid_spikes'], 'rt') as f:
    spikes = list(Bio.SeqIO.parse(f, 'fasta'))   
print(f"Read {len(spikes)} spike sequences.")
```

    Reading GISAID spikes in data/spikeprot0511.tar.tar.xz
    Read 1444763 spike sequences.


Make a data frame that has the BioPython SeqRecord, length, host, and geographic location (country) for each spike.
Also determine whether sequences have ambiguous amino acids or are all valid amino acids:


```python
spikes_df = (
    pd.DataFrame({'seqrecord': spikes})
    .assign(description=lambda x: x['seqrecord'].map(lambda rec: rec.description),
            country=lambda x: x['description'].str.split('|').str[-1],
            host=lambda x: x['description'].str.split('|').str[6].str.strip(),
            length=lambda x: x['seqrecord'].map(len),
            n_ambiguous=lambda x: x['seqrecord'].map(lambda rec: rec.seq.count('X') + rec.seq.count('x')),
            )
    )
```

Show number of sequences from different hosts, then keep only human ones:


```python
display(HTML(
    spikes_df
    .groupby('host')
    .aggregate(n_sequences=pd.NamedAgg('seqrecord', 'count'))
    .sort_values('n_sequences', ascending=False)
    .to_html()
    ))

spikes_df = spikes_df.query('host == "Human"')
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>n_sequences</th>
    </tr>
    <tr>
      <th>host</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Human</th>
      <td>1443052</td>
    </tr>
    <tr>
      <th>Neovison vison</th>
      <td>921</td>
    </tr>
    <tr>
      <th>Environment</th>
      <td>666</td>
    </tr>
    <tr>
      <th>Felis catus</th>
      <td>30</td>
    </tr>
    <tr>
      <th>Panthera leo</th>
      <td>24</td>
    </tr>
    <tr>
      <th>Manis javanica</th>
      <td>19</td>
    </tr>
    <tr>
      <th>Canis lupus familiaris</th>
      <td>18</td>
    </tr>
    <tr>
      <th>Panthera tigris jacksoni</th>
      <td>7</td>
    </tr>
    <tr>
      <th>Rhinolophus malayanus</th>
      <td>4</td>
    </tr>
    <tr>
      <th>Gorilla gorilla gorilla</th>
      <td>3</td>
    </tr>
    <tr>
      <th>Mus musculus</th>
      <td>3</td>
    </tr>
    <tr>
      <th>Rhinolophus shameli</th>
      <td>2</td>
    </tr>
    <tr>
      <th>Rhinolophus stheno</th>
      <td>2</td>
    </tr>
    <tr>
      <th>Mesocricetus auratus</th>
      <td>2</td>
    </tr>
    <tr>
      <th>Chlorocebus sabaeus</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Mus musculus (BALB/c mice)</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Prionailurus bengalensis euptilurus</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Rhinolophus affinis</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Rhinolophus bat</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Manis pentadactyla</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Rhinolophus pusillus</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Host</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Rhinolophus sinicus</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Mustela putorius furo</th>
      <td>1</td>
    </tr>
  </tbody>
</table>


Plot distribution of lengths and only keep sequences that are full-length (or near full-length) spikes:


```python
print('Distribution of length for all spikes:')
p = (ggplot(spikes_df) +
     aes('length') +
     geom_bar() +
     ylab('number of sequences') +
     theme(figure_size=(10, 2))
     )
fig = p.draw()
display(fig)
plt.close(fig)

min_length, max_length = 1260, 1276
print(f"\nOnly keeping spikes with lengths between {min_length} and {max_length}")
spikes_df = (
    spikes_df
    .assign(valid_length=lambda x: (min_length <= x['length']) & (x['length'] <= max_length))
    )

print('Here are number of sequences with valid and invalid lengths:')
display(HTML(spikes_df
             .groupby('valid_length')
             .aggregate(n_sequences=pd.NamedAgg('seqrecord', 'count'))
             .to_html()
             ))

print('\nDistribution of lengths for sequences with valid and invalid lengths; '
      'dotted red lines delimit valid lengths:')
p = (ggplot(spikes_df
            .assign(valid_length=lambda x: x['valid_length'].map({True: 'valid length',
                                                                  False: 'invalid length'}))
            ) +
     aes('length') +
     geom_bar() +
     ylab('number of sequences') +
     theme(figure_size=(10, 2), subplots_adjust={'wspace': 0.2}) +
     facet_wrap('~ valid_length', scales='free') +
     geom_vline(xintercept=min_length - 0.5, color='red', linetype='dotted') +
     geom_vline(xintercept=max_length + 0.5, color='red', linetype='dotted')
     )
fig = p.draw()
display(fig)
plt.close(fig)

spikes_df = spikes_df.query('valid_length')
```

    Distribution of length for all spikes:



    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_14_1.png)
    


    
    Only keeping spikes with lengths between 1260 and 1276
    Here are number of sequences with valid and invalid lengths:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>n_sequences</th>
    </tr>
    <tr>
      <th>valid_length</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>False</th>
      <td>28645</td>
    </tr>
    <tr>
      <th>True</th>
      <td>1414407</td>
    </tr>
  </tbody>
</table>


    
    Distribution of lengths for sequences with valid and invalid lengths; dotted red lines delimit valid lengths:



    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_14_5.png)
    


Finally, we get rid of spikes with **lots** of ambiguous residues as they may hinder the alignment below.
We will then do more detailed filtering for ambiguous residues just in the RBD region after alignment:


```python
max_ambiguous = 100
print(f"Filtering sequences with > {max_ambiguous} ambiguous residues")
spikes_df = (
    spikes_df
    .assign(excess_ambiguous=lambda x: x['n_ambiguous'] > max_ambiguous)
    )
display(HTML(
    spikes_df
    .groupby('excess_ambiguous')
    .aggregate(n_sequences=pd.NamedAgg('seqrecord', 'count'))
    .to_html()
    ))
```

    Filtering sequences with > 100 ambiguous residues



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>n_sequences</th>
    </tr>
    <tr>
      <th>excess_ambiguous</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>False</th>
      <td>1342655</td>
    </tr>
    <tr>
      <th>True</th>
      <td>71752</td>
    </tr>
  </tbody>
</table>


## Align the RBD region of the spikes
We now align the RBD regions of the spikes.
We do this **before** we filter sequences with too many ambiguous residues so that we can do that filtering just on the RBD region.

We align with `mafft` using the `--addfragments` and `--keeplength` options (see [here](https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html) and [here](https://mafft.cbrc.jp/alignment/software/addsequences.html)) to align relative to a reference that is just the RBD; these options also clip the sequences relative to the reference.
Note that these options make sense if the following conditions are met:
  1. Sequences are all very similar.
  2. We are not worried about insertions.
For now, both of these appear to be true, but this choice should be kept in mind if there is a lot more divergence or insertions.

We align relative to the reference that is the wildtype sequence used for the experiments:


```python
print(f"Reading reference nucleotide sequence in {config['wildtype_sequence']}")
refseq = Bio.SeqIO.read(config['wildtype_sequence'], 'fasta')

refprotfile = os.path.join(config['gisaid_mutations_dir'], 'reference_RBD.fasta')
print(f"Writing protein translation of reference sequence to {refprotfile}")
refseq.seq = refseq.seq.translate()
_ = Bio.SeqIO.write(refseq, refprotfile, 'fasta')
```

    Reading reference nucleotide sequence in data/wildtype_sequence.fasta
    Writing protein translation of reference sequence to results/GISAID_mutations/reference_RBD.fasta


Write all the other spikes to a file:


```python
spikes_file = os.path.join(config['gisaid_mutations_dir'],
                           'human_full-length_spikes.fasta')
print(f"Writing the spikes to {spikes_file}")
_ = Bio.SeqIO.write(spikes_df['seqrecord'].tolist(), spikes_file, 'fasta')
```

    Writing the spikes to results/GISAID_mutations/human_full-length_spikes.fasta


Now make the alignment.
Note that we use multiple threads to speed things up, and also align the spikes in chunks.
The reason that we have to the chunkwise alignment is that some unclear `mafft` error was arising if we tried to align them all at once:


```python
chunksize = 50000

aligned_rbds = []

for i in range(0, len(spikes_df), chunksize):
    spikes_file = os.path.join(config['gisaid_mutations_dir'],
                               f"human_full-length_spikes_{i + 1}-to-{i + chunksize}.fasta")
    print(f"Writing spikes {i + 1} to {i + chunksize} to {spikes_file}")
    _ = Bio.SeqIO.write(spikes_df['seqrecord'].tolist()[i: i + chunksize], spikes_file, 'fasta')
    print('Now aligning these sequences...')
    cmds = ['mafft', '--auto', '--thread', str(config['max_cpus']),
            '--keeplength', '--addfragments', spikes_file, refprotfile]
    res = subprocess.run(cmds, capture_output=True)
    if res.returncode:
        raise RuntimeError(f"Error in alignment:\n{res.stderr}")
    else:
        print('Alignment complete.\n')
        with io.StringIO(res.stdout.decode('utf-8')) as f:
            iseqs = list(Bio.SeqIO.parse(f, 'fasta'))
            # remove reference sequence, which should be first in file
            assert iseqs[0].seq == refseq.seq and iseqs[0].description == refseq.description
            iseqs = iseqs[1:]
            assert len(iseqs) == min(chunksize, len(spikes_df) - i)
            aligned_rbds += iseqs
            
assert len(aligned_rbds) == len(spikes_df)
```

    Writing spikes 1 to 50000 to results/GISAID_mutations/human_full-length_spikes_1-to-50000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 50001 to 100000 to results/GISAID_mutations/human_full-length_spikes_50001-to-100000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 100001 to 150000 to results/GISAID_mutations/human_full-length_spikes_100001-to-150000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 150001 to 200000 to results/GISAID_mutations/human_full-length_spikes_150001-to-200000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 200001 to 250000 to results/GISAID_mutations/human_full-length_spikes_200001-to-250000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 250001 to 300000 to results/GISAID_mutations/human_full-length_spikes_250001-to-300000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 300001 to 350000 to results/GISAID_mutations/human_full-length_spikes_300001-to-350000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 350001 to 400000 to results/GISAID_mutations/human_full-length_spikes_350001-to-400000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 400001 to 450000 to results/GISAID_mutations/human_full-length_spikes_400001-to-450000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 450001 to 500000 to results/GISAID_mutations/human_full-length_spikes_450001-to-500000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 500001 to 550000 to results/GISAID_mutations/human_full-length_spikes_500001-to-550000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 550001 to 600000 to results/GISAID_mutations/human_full-length_spikes_550001-to-600000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 600001 to 650000 to results/GISAID_mutations/human_full-length_spikes_600001-to-650000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 650001 to 700000 to results/GISAID_mutations/human_full-length_spikes_650001-to-700000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 700001 to 750000 to results/GISAID_mutations/human_full-length_spikes_700001-to-750000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 750001 to 800000 to results/GISAID_mutations/human_full-length_spikes_750001-to-800000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 800001 to 850000 to results/GISAID_mutations/human_full-length_spikes_800001-to-850000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 850001 to 900000 to results/GISAID_mutations/human_full-length_spikes_850001-to-900000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 900001 to 950000 to results/GISAID_mutations/human_full-length_spikes_900001-to-950000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 950001 to 1000000 to results/GISAID_mutations/human_full-length_spikes_950001-to-1000000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1000001 to 1050000 to results/GISAID_mutations/human_full-length_spikes_1000001-to-1050000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1050001 to 1100000 to results/GISAID_mutations/human_full-length_spikes_1050001-to-1100000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1100001 to 1150000 to results/GISAID_mutations/human_full-length_spikes_1100001-to-1150000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1150001 to 1200000 to results/GISAID_mutations/human_full-length_spikes_1150001-to-1200000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1200001 to 1250000 to results/GISAID_mutations/human_full-length_spikes_1200001-to-1250000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1250001 to 1300000 to results/GISAID_mutations/human_full-length_spikes_1250001-to-1300000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1300001 to 1350000 to results/GISAID_mutations/human_full-length_spikes_1300001-to-1350000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1350001 to 1400000 to results/GISAID_mutations/human_full-length_spikes_1350001-to-1400000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1400001 to 1450000 to results/GISAID_mutations/human_full-length_spikes_1400001-to-1450000.fasta
    Now aligning these sequences...
    Alignment complete.
    


## Parse / filter aligned RBDs

Now put all of the aligned RBDs into a data frame to filter and parse:


```python
rbd_df = (
    pd.DataFrame({'seqrecord': aligned_rbds})
    .assign(description=lambda x: x['seqrecord'].map(lambda rec: rec.description),
            country=lambda x: x['description'].str.split('|').str[-1],
            host=lambda x: x['description'].str.split('|').str[6].str.strip(),
            length=lambda x: x['seqrecord'].map(len),
            n_ambiguous=lambda x: x['seqrecord'].map(lambda rec: rec.seq.count('X') + rec.seq.count('x')),
            n_gaps=lambda x: x['seqrecord'].map(lambda rec: rec.seq.count('-')),
            all_valid_aas=lambda x: x['seqrecord'].map(lambda rec: re.fullmatch(f"[{protein_letters}]+",
                                                                                str(rec.seq)) is not None),
            )
    )

assert all(rbd_df['length'] == len(refseq))
```

Plot number of gaps and ambiguous nucleotides among sequences:


```python
for prop in ['n_ambiguous', 'n_gaps']:
    p = (ggplot(rbd_df) +
         aes(prop) +
         ylab('number of sequences') +
         theme(figure_size=(10, 2.5)) +
         geom_bar()
         )
    _ = p.draw()
```


    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_26_0.png)
    



    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_26_1.png)
    


Based on above plots, we will retain just RBDs with no ambiguous amino acids and no gaps:


```python
rbd_df = rbd_df.query('n_ambiguous == 0').query('n_gaps == 0')
assert rbd_df['all_valid_aas'].all()
print(f"Retained {len(rbd_df)} RBDs.")
```

    Retained 1327933 RBDs.


Now get and plot the number of amino-acid mutations per RBD relative to the reference sequence, plotting on both a linear and log scale.
We then filter all RBDs that have more than some maximum number of mutations, based on the idea that ones that are extremely highly mutated probably are erroneous.
**Note that this maximum number of mutations will change over time, so should be re-assessed periodically by looking at below plots.**


```python
max_muts = 8

refseq_str = str(refseq.seq)
rbd_df = (
    rbd_df
    .assign(seq=lambda x: x['seqrecord'].map(lambda rec: str(rec.seq)),
            n_mutations=lambda x: x['seq'].map(lambda s: sum(x != y for x, y in zip(s, refseq_str))))
    )

p = (ggplot(rbd_df) +
     aes('n_mutations') +
     geom_bar() +
     theme(figure_size=(10, 2.5)) +
     geom_vline(xintercept=max_muts + 0.5, color='red', linetype='dotted')
     )
_ = p.draw()
_ = (p + scale_y_log10()).draw()

rbd_df = rbd_df.query('n_mutations <= @max_muts')
```


    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_30_0.png)
    



    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_30_1.png)
    


Write RBD sequences that pass filtering to a file:


```python
print(f"Overall, there are {len(rbd_df)} aligned RBDs that passed filters.")

rbd_alignment_file = os.path.join(config['gisaid_mutations_dir'], 'RBD_alignment.fasta')
print(f"Writing alignment to {rbd_alignment_file}")
_ = Bio.SeqIO.write(rbd_df['seqrecord'].tolist(), rbd_alignment_file, 'fasta')
```

    Overall, there are 1327813 aligned RBDs that passed filters.
    Writing alignment to results/GISAID_mutations/RBD_alignment.fasta


## Get counts of each mutation
Now we get a data frame that gives the count of each mutation at each site:


```python
records = []
for tup in rbd_df[['seq', 'country']].itertuples():
    for isite, (mut, wt) in enumerate(zip(tup.seq, refseq_str), start=1):
        if mut != wt:
            records.append((isite, isite + config['site_number_offset'], wt, mut, tup.country))
            
muts_df = (pd.DataFrame.from_records(records,
                                     columns=['isite', 'site', 'wildtype', 'mutant', 'country'])
           .groupby(['isite', 'site', 'wildtype', 'mutant'])
           .aggregate(count=pd.NamedAgg('country', 'count'),
                      n_countries=pd.NamedAgg('country', 'nunique'))
           .reset_index()
           .sort_values('count', ascending=False)
           .assign(frequency=lambda x: x['count'] / len(rbd_df))
           )

print('Here are first few lines of mutation counts data frame:')
display(HTML(muts_df.head(n=15).to_html(index=False)))
```

    Here are first few lines of mutation counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>isite</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutant</th>
      <th>count</th>
      <th>n_countries</th>
      <th>frequency</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>171</td>
      <td>501</td>
      <td>N</td>
      <td>Y</td>
      <td>560239</td>
      <td>470</td>
      <td>0.421926</td>
    </tr>
    <tr>
      <td>154</td>
      <td>484</td>
      <td>E</td>
      <td>K</td>
      <td>50055</td>
      <td>253</td>
      <td>0.037697</td>
    </tr>
    <tr>
      <td>122</td>
      <td>452</td>
      <td>L</td>
      <td>R</td>
      <td>49525</td>
      <td>130</td>
      <td>0.037298</td>
    </tr>
    <tr>
      <td>147</td>
      <td>477</td>
      <td>S</td>
      <td>N</td>
      <td>46922</td>
      <td>181</td>
      <td>0.035338</td>
    </tr>
    <tr>
      <td>109</td>
      <td>439</td>
      <td>N</td>
      <td>K</td>
      <td>19861</td>
      <td>163</td>
      <td>0.014958</td>
    </tr>
    <tr>
      <td>148</td>
      <td>478</td>
      <td>T</td>
      <td>K</td>
      <td>15175</td>
      <td>72</td>
      <td>0.011429</td>
    </tr>
    <tr>
      <td>87</td>
      <td>417</td>
      <td>K</td>
      <td>N</td>
      <td>12354</td>
      <td>137</td>
      <td>0.009304</td>
    </tr>
    <tr>
      <td>87</td>
      <td>417</td>
      <td>K</td>
      <td>T</td>
      <td>10070</td>
      <td>83</td>
      <td>0.007584</td>
    </tr>
    <tr>
      <td>164</td>
      <td>494</td>
      <td>S</td>
      <td>P</td>
      <td>6469</td>
      <td>76</td>
      <td>0.004872</td>
    </tr>
    <tr>
      <td>190</td>
      <td>520</td>
      <td>A</td>
      <td>S</td>
      <td>3320</td>
      <td>68</td>
      <td>0.002500</td>
    </tr>
    <tr>
      <td>171</td>
      <td>501</td>
      <td>N</td>
      <td>T</td>
      <td>2938</td>
      <td>71</td>
      <td>0.002213</td>
    </tr>
    <tr>
      <td>192</td>
      <td>522</td>
      <td>A</td>
      <td>S</td>
      <td>2207</td>
      <td>79</td>
      <td>0.001662</td>
    </tr>
    <tr>
      <td>154</td>
      <td>484</td>
      <td>E</td>
      <td>Q</td>
      <td>1978</td>
      <td>52</td>
      <td>0.001490</td>
    </tr>
    <tr>
      <td>160</td>
      <td>490</td>
      <td>F</td>
      <td>S</td>
      <td>1892</td>
      <td>75</td>
      <td>0.001425</td>
    </tr>
    <tr>
      <td>110</td>
      <td>440</td>
      <td>N</td>
      <td>K</td>
      <td>1783</td>
      <td>62</td>
      <td>0.001343</td>
    </tr>
  </tbody>
</table>


Plot how many mutations are observed how many times:


```python
p = (ggplot(muts_df) +
     aes('count') +
     geom_histogram(bins=20) +
     scale_x_log10() +
     ylab('number of sequences') +
     xlab('times mutation observed')
     )

_ = p.draw()
```


    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_36_0.png)
    


Write the mutation counts to a file:


```python
print(f"Writing mutation counts to {config['gisaid_mutations_dir']}")
muts_df.to_csv(config['gisaid_mutation_counts'], index=False)
```

    Writing mutation counts to results/GISAID_mutations

