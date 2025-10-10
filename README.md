# Nanopore AMR and Pangenome Analysis Pipeline

A comprehensive Nextflow pipeline for bacterial genome analysis from Nanopore sequencing data or pre-assembled genomes. This pipeline performs quality control, assembly, antimicrobial resistance (AMR) gene detection, plasmid identification, MLST typing, functional annotation, and pangenome analysis.

## Features

- **Quality Control**: NanoPlot QC reports for raw and filtered reads
- **Read Filtering**: Filtlong for quality and length filtering
- **Assembly**: Flye assembler optimized for Nanopore reads
- **Assembly Evaluation**: QUAST metrics and Bandage visualization
- **AMR Detection**: ABRicate screening against multiple databases (ResFinder, PlasmidFinder)
- **MLST Typing**: Multi-locus sequence typing for strain identification
- **Genome Annotation**: Bakta for comprehensive functional annotation
- **Pangenome Analysis**: Panaroo for comparative genomics
- **Phylogenetics**: Optional RAxML-NG and IQ-TREE phylogenetic trees
- **Taxonomy**: Kraken2 and GTDB-Tk taxonomic classification
- **Reporting**: MultiQC aggregated reports

## Quick Start

### Prerequisites

- Nextflow (≥ 21.04)
- Conda/Mamba (for software dependencies)
- Reference databases:
  - Bakta database (required for annotation)
  - Kraken2 database (optional, for taxonomy)
  - GTDB-Tk database (optional, for taxonomy)

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/nanopore-amr-pangenome.git
cd nanopore-amr-pangenome

# Test the pipeline
nextflow run main.nf --help
```

### Basic Usage

```bash
# Run with Nanopore reads (full pipeline with pangenome)
nextflow run main.nf --reads 'data/*.fastq.gz' --bakta_db /path/to/bakta_db

# Run with pre-assembled genomes
nextflow run main.nf -entry amr_annotation_from_assemblies \
  --assemblies 'assemblies/*.fasta' \
  --bakta_db /path/to/bakta_db
```

## Pipeline Workflows

The pipeline offers five workflows optimized for different use cases:

### 1. `assembly_amr_pangenome` (default)

**Complete pipeline from raw reads to pangenome analysis**

**Input:** Nanopore FASTQ files  
**Output:** QC, assemblies, AMR/plasmids, MLST, annotation, pangenome

**Steps:**
1. Quality control (NanoPlot)
2. Read filtering (Filtlong)
3. Assembly (Flye)
4. Assembly evaluation (QUAST, Bandage)
5. AMR/plasmid detection (ABRicate)
6. MLST typing
7. Genome annotation (Bakta)
8. Pangenome analysis (Panaroo)
9. GML visualization

```bash
nextflow run main.nf \
  --reads 'data/*.fastq.gz' \
  --bakta_db /path/to/bakta_db
```

**Use when:** You have raw Nanopore reads and want comprehensive analysis including pangenome

---

### 2. `assembly_amr_annotation`

**Individual sample analysis without pangenome**

**Input:** Nanopore FASTQ files  
**Output:** QC, assemblies, AMR/plasmids, MLST, annotation (no pangenome)

**Steps:** Same as workflow 1, but stops after Bakta annotation

```bash
nextflow run main.nf -entry assembly_amr_annotation \
  --reads 'data/*.fastq.gz' \
  --bakta_db /path/to/bakta_db
```

**Use when:** 
- You need rapid screening without pangenome overhead
- Analyzing samples independently
- Resource-constrained environments
- Preliminary quality assessment

---

### 3. `amr_pangenome_from_assemblies`

**Analysis from pre-assembled genomes with pangenome**

**Input:** Pre-assembled FASTA files  
**Output:** Assembly evaluation, AMR/plasmids, MLST, annotation, pangenome

**Steps:**
1. Assembly evaluation (QUAST)
2. AMR/plasmid detection (ABRicate)
3. MLST typing
4. Genome annotation (Bakta)
5. Pangenome analysis (Panaroo)
6. GML visualization

```bash
nextflow run main.nf -entry amr_pangenome_from_assemblies \
  --assemblies 'assemblies/*.fasta' \
  --bakta_db /path/to/bakta_db
```

**Use when:**
- You have pre-assembled genomes from any source
- Working with external/public assemblies
- Need comparative genomics on existing assemblies

---

### 4. `amr_annotation_from_assemblies`

**Individual assembly analysis without pangenome**

**Input:** Pre-assembled FASTA files  
**Output:** Assembly evaluation, AMR/plasmids, MLST, annotation (no pangenome)

**Steps:** Same as workflow 3, but stops after Bakta annotation

```bash
nextflow run main.nf -entry amr_annotation_from_assemblies \
  --assemblies 'assemblies/*.fasta' \
  --bakta_db /path/to/bakta_db
```

**Use when:**
- Screening large numbers of assemblies (100+)
- Quick AMR/plasmid profiling
- Pre-checking assembly quality
- Mixed-source assemblies from different methods

---

### 5. `taxonomy_from_reads`

**Taxonomic classification from reads**

**Input:** Nanopore FASTQ files  
**Output:** Kraken2 classification, assembly, GTDB-Tk taxonomy

**Steps:**
1. Kraken2 read classification
2. Assembly (Flye)
3. GTDB-Tk genome classification

```bash
nextflow run main.nf -entry taxonomy_from_reads \
  --reads 'data/*.fastq.gz' \
  --kraken2_db /path/to/kraken2_db \
  --gtdbtk_db /path/to/gtdbtk_db
```

**Use when:** 
- Confirming species identification
- Checking for contamination
- Unknown or mixed samples

---

## Workflow Decision Tree

```
Do you have raw reads or assemblies?
│
├─ Raw Reads (FASTQ)
│  │
│  └─ Do you need pangenome analysis?
│     ├─ YES → assembly_amr_pangenome (default)
│     ├─ NO  → assembly_amr_annotation
│     └─ Just taxonomy → taxonomy_from_reads
│
└─ Assemblies (FASTA)
   │
   └─ Do you need pangenome analysis?
      ├─ YES → amr_pangenome_from_assemblies
      └─ NO  → amr_annotation_from_assemblies
```

## Configuration

### Essential Parameters

Edit `nextflow.config` or provide via command line:

```bash
# Input/Output
--reads          'path/to/*.fastq.gz'    # Nanopore reads
--assemblies     'path/to/*.fasta'       # Pre-assembled genomes
--outdir         'results'               # Output directory

# Databases (required)
--bakta_db       '/path/to/bakta/db'     # Bakta annotation database

# ABRicate databases
--abricate_db_amr  'resfinder'           # AMR database
--abricate_db_plm  'plasmidfinder'       # Plasmid database

# MLST
--mlst_scheme    'ecoli_achtman_4'       # MLST scheme

# Taxonomy (optional)
--kraken2_db     '/path/to/kraken2/db'   # Kraken2 database
--gtdbtk_db      '/path/to/gtdbtk/db'    # GTDB-Tk database

# Panaroo
--panaroo_clean_mode  'strict'           # strict, moderate, sensitive
--panaroo_plot_max_nodes  800            # Max nodes to plot
```

### Resource Profiles

The pipeline includes resource profiles for different environments:

```bash
# Local execution (default)
nextflow run main.nf --reads 'data/*.fastq.gz'

# SLURM cluster
nextflow run main.nf -profile slurm --reads 'data/*.fastq.gz'

# With Docker
nextflow run main.nf -profile docker --reads 'data/*.fastq.gz'

# With Singularity
nextflow run main.nf -profile singularity --reads 'data/*.fastq.gz'
```

### Process Labels

Processes are labeled by resource requirements:

- **light**: 2 CPUs, 4 GB RAM, 2 hours (QC, filtering, small tools)
- **medium**: 8 CPUs, 16 GB RAM, 8 hours (assembly, annotation)
- **heavy**: 16 CPUs, 48 GB RAM, 24 hours (pangenome, phylogenetics)

## Output Structure

```
results/
├── nanoplot_raw/           # Raw read QC reports
├── filtlong/               # Filtered reads
├── nanoplot_filt/          # Filtered read QC reports
├── flye/                   # Assemblies and assembly graphs
├── quast/                  # Assembly quality metrics
├── bandage/                # Assembly graph visualizations
├── abricate/               # AMR and plasmid gene results
│   ├── *_amr.tsv          # Per-sample AMR genes
│   ├── *_plm.tsv          # Per-sample plasmid genes
│   └── abricate_merged_all.tsv  # Combined results
├── abricate_summary/       # AMR/plasmid summary plots
├── mlst/                   # MLST typing results
│   ├── *.mlst.tsv         # Per-sample MLST
│   ├── mlst_merged_all.tsv
│   └── mlst_summary_counts.tsv
├── bakta/                  # Genome annotations
│   ├── *.gff3             # Gene features
│   ├── *.faa              # Protein sequences
│   ├── *.gbff             # GenBank format
│   └── *.png              # Genome plots
├── panaroo/                # Pangenome analysis (if enabled)
│   ├── results/
│   │   ├── core_gene_alignment.aln.fasta
│   │   ├── gene_presence_absence.csv
│   │   └── final_graph.gml
│   └── plots/             # Network visualizations
├── multiqc/                # Aggregated QC report
├── kraken2/                # Taxonomic classification (if run)
├── gtdbtk/                 # GTDB taxonomy (if run)
├── pipeline_report.html    # Nextflow execution report
├── pipeline_timeline.html  # Execution timeline
└── pipeline_dag.png        # Workflow diagram
```

## Key Output Files

| File | Description |
|------|-------------|
| `multiqc/multiqc_report.html` | Aggregated QC and results summary |
| `abricate/abricate_merged_all.tsv` | All AMR and plasmid genes detected |
| `abricate_summary/*.png` | Visualizations of AMR/plasmid profiles |
| `mlst/mlst_summary_counts.tsv` | MLST sequence type distribution |
| `bakta/*/*.gff3` | Genome annotations in GFF3 format |
| `panaroo/results/gene_presence_absence.csv` | Pangenome gene matrix |
| `panaroo/plots/LCC_plot.png` | Pangenome network visualization |

## Feature Comparison

| Feature | With Pangenome | Without Pangenome |
|---------|----------------|-------------------|
| NanoPlot QC | ✅ | ✅ |
| Filtlong | ✅ | ✅ |
| Flye Assembly | ✅ | ✅ |
| QUAST Evaluation | ✅ | ✅ |
| Bandage Visualization | ✅ | ✅ |
| ABRicate (AMR/Plasmid) | ✅ | ✅ |
| MLST Typing | ✅ | ✅ |
| Bakta Annotation | ✅ | ✅ |
| MultiQC Report | ✅ | ✅ |
| **Core Gene Alignment** | ✅ | ❌ |
| **Gene Presence/Absence** | ✅ | ❌ |
| **Pangenome Graph** | ✅ | ❌ |
| **GML Visualizations** | ✅ | ❌ |
| **Phylogenetic Trees** | ✅ | ❌ |

## Database Setup

### Bakta Database (Required)

```bash
# Download Bakta database
bakta_db download --output /path/to/bakta_db --type light
```

### Kraken2 Database (Optional)

```bash
# Download standard database
kraken2-build --standard --db /path/to/kraken2_db

# Or use pre-built databases
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220926.tar.gz
tar -xvzf k2_standard_20220926.tar.gz -C /path/to/kraken2_db
```

### GTDB-Tk Database (Optional)

```bash
# Download GTDB-Tk reference data
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar -xvzf gtdbtk_data.tar.gz -C /path/to/gtdbtk_db
```

## Advanced Usage

### Custom ABRicate Databases

```bash
# List available databases
abricate --list

# Use custom database
nextflow run main.nf \
  --reads 'data/*.fastq.gz' \
  --abricate_db_amr 'ncbi' \
  --abricate_db_plm 'plasmidfinder'
```

### Resume Failed Runs

Nextflow caches completed tasks, allowing you to resume:

```bash
nextflow run main.nf --reads 'data/*.fastq.gz' -resume
```

### Clean Work Directory

```bash
# Remove work directory after successful completion
nextflow run main.nf --reads 'data/*.fastq.gz'
nextflow clean -f
```

### Generate Only Workflow Diagram

```bash
nextflow run main.nf --reads 'test.fastq.gz' -with-dag flowchart.png -preview
```

## Troubleshooting

### Common Issues

**Issue**: `No reads found for: data/*.fastq.gz`  
**Solution**: Check file paths and use absolute paths or proper wildcards

**Issue**: Bakta fails with "Database not found"  
**Solution**: Verify `--bakta_db` path points to uncompressed database directory

**Issue**: Panaroo fails with "Invalid GFF"  
**Solution**: Ensure Bakta completed successfully; check `.gff3` files are not empty

**Issue**: MLST returns no results  
**Solution**: Verify `--mlst_scheme` matches your species; use `mlst --list` to see schemes

**Issue**: Out of memory errors  
**Solution**: Adjust resources in `nextflow.config` or use `-profile slurm` for HPC

### Getting Help

```bash
# View pipeline help
nextflow run main.nf --help

# Check specific workflow
nextflow run main.nf -entry assembly_amr_annotation --help

# Validate configuration
nextflow config main.nf
```

## Performance Tips

1. **Use `-resume`** to restart from cached results after failures
2. **Increase `maxForks`** in config for parallel execution on HPC
3. **Skip pangenome** for large datasets (100+ samples) to save resources
4. **Use `--bakta_skip_crispr`** on ARM/Mac to avoid PilerCR issues
5. **Pre-filter reads** if you have high-depth data (>100X coverage)

## Citation

If you use this pipeline, please cite the tools it uses:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35, 316–319.
- **Flye**: Kolmogorov, M., et al. (2019). Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology, 37, 540–546.
- **ABRicate**: Seemann T. ABRicate. https://github.com/tseemann/abricate
- **Bakta**: Schwengers, O., et al. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics, 7(11).
- **Panaroo**: Tonkin-Hill, G., et al. (2020). Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biology, 21, 180.
- **QUAST**: Gurevich, A., et al. (2013). QUAST: quality assessment tool for genome assemblies. Bioinformatics, 29(8), 1072-1075.
- **NanoPlot**: De Coster, W., et al. (2018). NanoPack: visualizing and processing long-read sequencing data. Bioinformatics, 34(15), 2666-2669.

## License

This pipeline is distributed under the MIT License. See `LICENSE` file for details.

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with clear description

## Contact

For questions, issues, or suggestions:
- Open an issue on GitHub
- Email: your.email@institution.edu

## Changelog

### Version 1.0.0
- Initial release
- Five workflow options for different use cases
- Support for both reads and assemblies
- Comprehensive AMR, MLST, and annotation
- Optional pangenome analysis with visualization
