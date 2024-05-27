### GTDB Strain downloader

Providing one or multiple species name, this pipeline will:
* Retrieve from NCBI all strain genomes related to this species based on GTDB metadata
* Run [SKANI](https://github.com/bluenote-1577/skani) and cluster genome at a strain level (ANI threshold = 1).

#### Requirements

Everything can be downloaded with conda:
* SKANI `conda install -c bioconda skani`
* SCIKIT-LEARN `conda install -c conda-forge scikit-learn`
* NCBI Datasets tool `conda install -c conda-forge ncbi-datasets-cli`

#### Configuration

You can edit the GTDB metadata file, the conda environments associated to each software or library (default `gtdb_sdl`), and the species of interest directly at the begining of the `gtdb_sdl.snk` file. For the `SPECIES` dictionnary, keys represent system names that will be reflected in the path, and values the actual name that will be searched in GTDB metadata.

I recommend getting a NCBI API key as well, but this is not mandatory.

The pipeline can then be run this way (`-n` is for dry-run, you can remove it):
```
snakemake -s gtdb_sdl.snk -c {threads} --use-conda -p --rerun-triggers mtime -n
```

The pipeline should take less than 10 minutes even with just 2 cores for each species.

#### OUTPUT

Outfiles can be found in the `output` folder where the snakemake pipeline has been launched. For each species, it contains:
* `{species}.genomes.tar` that contains genome for all reference strain.
* `{species}.meta.tsv.gz` that contains the metadata for all strains.
* `{species}.skani.tsv.gz` that contains the pairwise skani results.