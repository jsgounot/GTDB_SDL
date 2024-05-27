# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2024-05-27 11:23:48
# @Last Modified by:   jsgounot
# @Last Modified time: 2024-05-27 13:03:32

import pandas as pd

def extract_tax(values, target_rank):
    values = values.split(';')
    for value in values:
        rank, value = value.split('__')
        if rank == target_rank:
            return value
    return np.nan

fname = snakemake.input[0]
species = snakemake.params['sname']

df = []
chunksize = 10 ** 5

for chunk in pd.read_csv(fname, chunksize=chunksize, sep='\t'):
	chunk['species'] = chunk['gtdb_taxonomy'].apply(lambda tax: extract_tax(tax, 's'))
	chunk = chunk[chunk['species'] == species]
	if not chunk.empty: df.append(chunk)

if not df:
    raise Exception(f'No GTDB results were found with species name "{species}"')

df = pd.concat(df)

completeness = snakemake.params['completeness']
contamination = snakemake.params['contamination']
df = df[(df['checkm2_completeness'] >= completeness) & (df['checkm2_contamination'] <= contamination)]

outfile = snakemake.output['meta']
df.to_csv(outfile, sep='\t')

accessions = sorted(acc[3:] if (acc.startswith('GB_') or acc.startswith('RS_')) else acc 
    for acc in df['accession'])

with open(snakemake.output['acc'], 'w') as f:
	f.write('\n'.join(accessions))