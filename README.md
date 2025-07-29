# Nanopore AMR and Species Identification Pipeline

A modular, reproducible Nextflow DSL2 workflow for identifying bacterial species and AMR genes from Nanopore whole genome sequencing (WGS) data.

## Features
1. Modular DSL2 pipeline with process-specific scripts

2. Species classification via GTDB-Tk, Kraken2, FastANI

3. AMR gene detection using Abricate (ResFinder only)

4. Summary report with visual output (barplots, CSVs)

5. GitHub CI-tested with dummy data

6. Reproducible via environment.yml

## Quickstart
```
# Option 2: Clone and enter
git clone https://github.com/yourusername/nanopore_amr_pipeline.git
cd nanopore_amr_pipeline

# Options 2: Create conda env
mamba env create -f environment.yml
conda activate nanopore_amr

# Run the test dataset
nextflow run main.nf \
  --reads "data/*.fastq.gz" \
  --kraken_db /path/to/kraken_db \
  --gtdbtk_db /path/to/gtdbtk_db \
  --fastani_db /path/to/fastani_db \
  --outdir results
```
