# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-05-11 09:26:31
# @Last Modified by:   jsgounot
# @Last Modified time: 2024-05-27 16:01:01

import os, json
import pandas as pd
import numpy as np

fname = snakemake.input['mat']
df = pd.read_csv(fname, sep='\t')

# ----------------------------------------------------------------------
# Clustering
#   + scikit-learn     1.2.0  py310h6a678d5_1  anaconda/linux-64       9 MB

from sklearn.cluster import AgglomerativeClustering as AC
threshold = snakemake.params['ani_threshold']

df['ref_sample'] = df['Ref_file'].apply(os.path.basename)
df['que_sample'] = df['Query_file'].apply(os.path.basename)
df = df[['ref_sample', 'que_sample', 'ANI']]

dfr = df.copy()
dfr['ref_sample'] = df['que_sample']
dfr['que_sample'] = df['ref_sample']

dfexact = pd.DataFrame([
    {
        'ref_sample': name,
        'que_sample': name,
        'ANI': 100
    }
    for name in set(df['que_sample']) | set(df['ref_sample'])
])

df = pd.concat([df, dfr, dfexact])
matrix = df.pivot_table(index='ref_sample', columns='que_sample', values='ANI')
matrix = 100 - matrix

clustering = AC(linkage="average", n_clusters=None, compute_full_tree=True, 
                distance_threshold=threshold, metric="precomputed").fit_predict(matrix)

gidx = {genome: idx for idx, genome in enumerate(matrix.index)}
df['ref_cluster'] =  df['ref_sample'].map(gidx)
df['que_cluster'] =  df['que_sample'].map(gidx)

i2c = {idx: clusterID for idx, clusterID in enumerate(clustering)}
df['ref_cluster'] =  df['ref_cluster'].map(i2c)
df['que_cluster'] =  df['que_cluster'].map(i2c)

# ----------------------------------------------------------------------
# Cleaning and annotations

# Map path to accession number
fname = snakemake.input['ncbi_cat']
with open(fname) as f:
    jdata = json.load(f)

accession_to_filepath = atp = {
    line['accession']: [subline['filePath'] + '.gz' for subline in line['files']]
    for line in jdata['assemblies']
    if 'accession' in line
}

atpr = {}
for key, value in atp.items():
    assert len(value) == 1
    bname = os.path.basename(value[0])
    atpr[bname] = key

assert all(sample in atpr for sample in df['ref_sample'].unique())
assert all(sample in atpr for sample in df['que_sample'].unique())

df['ref_sample'] = df['ref_sample'].map(atpr)
df['que_sample'] = df['que_sample'].map(atpr)

outfile = snakemake.output['str_ani']
df.to_csv(outfile, sep='\t', index=False)

# Map metadata to results
fname = snakemake.input['meta']
mdata = pd.read_csv(fname, sep='\t', compression='gzip', index_col=0)

mdata = mdata[['accession', 'checkm2_completeness', 'checkm2_contamination', 'n50_scaffolds', 'gtdb_representative']]
mdata.columns = ['Genome', 'Completeness', 'Contamination', 'N50', 'Rep']
fun = lambda gid: gid[3:] if gid.startswith('RS_') or gid.startswith('GB_') else gid
mdata['Genome'] = mdata['Genome'].apply(fun)

# Merge and process
cdf = df[['ref_sample', 'ref_cluster']].drop_duplicates()
cdf.columns = ['Genome', 'ClusterID']
assert len(cdf) == len(set(df['ref_sample']) | set(df['que_sample']))
cdf = cdf.merge(mdata, on='Genome', how='left')

# Using UHGG score definition for Species representative
cdf['Score'] = cdf['Completeness'] - 5 * cdf['Contamination'] + 0.5 * np.log(cdf['N50'])

def select_rep(sdf):
    sdf = sdf.sort_values('Score', ascending=False)
    # We always favor the species representative if any
    rep = list(sdf[sdf['Rep'] == 't']['Genome'])
    rep = rep[0] if rep else list(sdf['Genome'])[0]
    sdf['Str_Rep'] = sdf['Genome'] == rep
    return sdf

cdf = cdf.groupby('ClusterID').apply(select_rep).reset_index(drop=True)

outfile = snakemake.output['str_meta']
cdf.to_csv(outfile, sep='\t', index=False)

reps = cdf[cdf['Str_Rep']]['Genome']
paths = [atp[rep][0] for rep in reps]

tpaths = [f'wdir/{snakemake.wildcards.species}/genomes/ncbi_dataset/data/{bname}'
    for bname in paths]

for path in tpaths:
    assert os.path.isfile(path)

outfile = snakemake.output['str_rlist']
with open(outfile, 'w') as f:
    f.write('\n'.join(paths))
