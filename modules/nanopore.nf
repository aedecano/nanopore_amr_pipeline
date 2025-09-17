#!/usr/bin/env nextflow

/*
#==============================================
Enable DSL2
#==============================================
*/

nextflow.enable.dsl = 2

/*
#==============================================
QC of raw reads
#==============================================
*/

process BASECALL_DORADO {

    label 'dorado basecall'
    tag {"Dorado basecalling ${sample_id} reads"}
    
    publishDir "$params.outdir/basecalled_dorado", mode: 'copy'

    label 'dorado basecall'
    tag   "${sample_id}"

    input:
    tuple val(sample_id), path(raw_dir)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    publishDir "${params.outdir}/DORADO_BASECALL", mode: 'copy', overwrite: true

    script:
    """
    dorado basecaller \\
        --device ${params.dorado_device} \\
        --emit-bam \\
        --emit-moves \\
        ${params.dorado_model} \\
        ${raw_dir} \\
        > ${sample_id}.bam
    """
}

process BAM_TO_FASTQ {
    label 'bam_to_fastq'
    tag   "${sample_id}"

    input:
    tuple val(sample_id), path(bam)

    output:
    path "${sample_id}.fastq.gz"

    publishDir "${params.outdir}/BAM_TO_FASTQ", mode: 'copy', overwrite: true

    script:
    """
    samtools fastq -n -F 0x900 ${bam} | gzip -c > ${sample_id}.fastq.gz
    """
}

process BASECALL_GUPPY {
    label 'guppy basecall'
    tag {"Guppy basecalling ${sample_id} reads"}
    
    publishDir "$params.outdir/basecalled", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) 
    
    output:
    tuple val(sample_id), path("${sample_id}_fastq/${sample_id}_1.fastq.gz"), path("${sample_id}_fastq/${sample_id}_2.fastq.gz")
    
    script:
    
    """
    guppy_basecaller -i ${reads} -s ${sample_id}_fastq --flowcell FLO-MIN106 --kit SQK-LSK109 --num_callers ${task.cpus} --cpu_threads_per_caller 4 --barcode_kits EXP-NBD104 --trim_barcodes
    """
}


process RAWFASTQC_SINGLE {
    label 'fastqc'
    tag {"FastQC raw ${sample_id} reads"}
    
    publishDir "$params.outdir/raw_fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) 
    
    output:
    //path ("${sample_id}_fastqc/${sample_id}.txt")
    path "${sample_id}_fastqc/*", emit: fastqc
    
    script:
    
    """
    fastqc --threads ${task.cpus} ${reads}
    """

}

process MULTIQC_READS {
    label 'multiqc'
    tag {"MultiQC raw reads"}
    
    publishDir "$params.outdir/multiqc_reads", mode: 'copy'

    input:
    path(reads) from RAWFASTQC_SINGLE.out.collect()
    
    output:
    path ("multiqc_report.html")
    
    script:
    
    """
    multiqc -o . ${reads}
    """

}

process FILTLONG {
  label 'filtlong'
  tag    "$sample_id"
  publishDir "${params.outdir}/filtlong", mode: 'copy'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}.filt.fastq.gz"), emit: reads

  // Use script: with Groovy interpolation; avoid `set -u`
  script:
  """
  set -eo pipefail
  echo "DEBUG sample_id=${sample_id}" >&2
  echo "DEBUG reads=${reads}" >&2
  test -s "${reads}" || { echo "Input reads missing/empty: ${reads}" >&2; exit 1; }

  filtlong --min_length 1000 --target_bases 200000000 "${reads}" \
    | gzip -c > "${sample_id}.filt.fastq.gz"
  """
}

process CLEANFASTQC_SINGLE {
    label 'fastqc'
    tag {"FastQC cleaned ${sample_id} reads"}
    
    publishDir "$params.outdir/cleaned_fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) 
    
    output:
    //path ("${sample_id}_clean_fastqc/${sample_id}.txt")
    path "${sample_id}_clean_fastqc/*", emit: fastqc
    
    script:
    
    """
    fastqc --threads ${task.cpus} ${reads}
    """

}      

process MULTIQC_CONTIGS {
    label 'multiqc'
    tag {"MultiQC contigs"}
    
    publishDir "$params.outdir/multiqc_contigs", mode: 'copy'
    input:
    path(reads) from CLEANFASTQC_SINGLE.out.collect()           
    output:
    path ("multiqc_report.html")
    script:
    
    """
    multiqc -o . ${reads}
    """ 
}

process FLYE {
  tag "!{sample_id}"
    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_assembly/${sample_id}_contigs.fasta"), emit: assembly
  
    script:      
        """
        flye --nano-raw "!{reads}" --out-dir . --threads ${task.cpus}
        mv assembly.fasta "!{sample_id}.contigs.fasta"
        """
}


process QUAST_FROM_READS {
    label 'quast'
    tag {"QUAST from ${sample_id} reads"}
    
    publishDir "$params.outdir/quast_reads", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) 
    
    output:
    path ("${sample_id}_quast_report")
    
    script:
    
    """
    quast.py ${reads} -o ${sample_id}_quast_report -t ${task.cpus}
    """
}   

process QUAST_FROM_CONTIGS {
    label 'quast'
    tag {"QUAST from ${sample_id} contigs"}
    
    publishDir "$params.outdir/quast_contigs", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs) 
    
    output:
    path ("${sample_id}_quast_report")
    
    script:
    
    """
    quast.py ${contigs} -o ${sample_id}_quast_report -t ${task.cpus}
    """
}   

process MEDAKA {
    label 'medaka'
    tag {"Medaka polishing ${sample_id} reads"}
    
    publishDir "$params.outdir/polished_contigs", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(contigs) 
    
    output:
    tuple val(sample_id), path("${sample_id}_polished/${sample_id}_polished.fasta")
    
    script:
    
    """
    medaka_consensus -i ${reads} -d ${contigs} -o ${sample_id}_polished -t ${task.cpus} -m r941_min_high_g360
    """
}

process NANOPLOTYPER {
    label 'nanoplot'
    tag {"NanoPlot ${sample_id} reads"}
    
    publishDir "$params.outdir/nanoplot", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) 
    
    output:
    path ("${sample_id}_nanoplot")
    
    script:
    
    """
    NanoPlot --fastq ${reads} -o ${sample_id}_nanoplot --threads ${task.cpus}
    """
}

process FASTANI {
    label 'fastani'
    tag {"FastANI ${sample_id} contigs"}
    
    publishDir "$params.outdir/fastani", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs) 
    
    output:
    path ("${sample_id}_fastani_report.txt")
    
    script:
    
    """
    fastANI -q ${contigs} -r ${ref}.fasta -o ${sample_id}_fastani_report.txt
    """
}

process GTDBTYPER {
    label 'gtdbtyper'
    tag {"GTDB-Tk classification ${sample_id} contigs"}
    
    publishDir "$params.outdir/gtdbtk", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs) 
    
    output:
    path ("${sample_id}_gtdbtk_report.txt")
    
    script:
    
    """
    gtdbtk classify_wf --genome_dir . --out_dir ${sample_id}_gtdbtk --cpus ${task.cpus}
    cp ${sample_id}_gtdbtk/gtdbtk.bac120.summary.tsv ${sample_id}_gtdbtk_report.txt
    """
}


process AMR_PLM_FROM_CONTIGS {
    label 'amr_plm'
    tag {"AMR Plasmid from ${sample_id} contigs"}
    
    publishDir "$params.outdir/amr_plasmids_contigs", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs) 
    
    output:
    path ("${sample_id}_plasmid_report.txt")
    
    script:
    
    """
    plasmidfinder.py -i ${contigs} -o ${sample_id}_plasmid_report -p /path/to/plasmidfinder_db -t 0.95 -l 0.6
    """
}

process AMR_ABRFORMAT {
  tag "${sample_id}"

  publishDir "results/amr", mode: 'copy', overwrite: true, pattern: "*.txt"

  input:
    tuple val(sample_id), path(contigs)

  output:
    path "${sample_id}_abricate_report.txt", emit: abricate_report

  script:
  """
  abricate --db resfinder ${contigs} > ${sample_id}_abricate_report.txt
  """
}

process MOBTYPER {
    label 'mobtyper'
    tag {"MOB-typer from ${sample_id} contigs"}
    
    publishDir "$params.outdir/mobtyper_contigs", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs) 
    
    output:
    path ("${sample_id}_mobtyper_report.txt")
    
    script:
    
    """
    mob_typer -i ${contigs} -o ${sample_id}_mobtyper_report
    """
}

process AMRFINDERPLUS {
    label 'amrfinderplus'
    tag {"AMR Finder Plus from ${sample_id} contigs"}
    
    publishDir "$params.outdir/amrfinderplus_contigs", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs) 
    
    output:
    path ("${sample_id}_amrfinderplus_report.txt")
    
    script:
    
    """
    amrfinder -n ${contigs} -o ${sample_id}_amrfinderplus_report.txt --organism Clostridioides_difficile --plus
    """
}

process MLST_FROM_CONTIGS {
    label 'mlst'
    tag {"MLST from ${sample_id} contigs"}
    
    publishDir "$params.outdir/mlst_contigs", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs) 
    
    output:
    path ("${sample_id}_mlst_report.txt")
    
    script:
    
    """
    mlst ${contigs} > ${sample_id}_mlst_report.txt
    """
}

process PROKKA {
    label 'prokka'
    tag {"Prokka annotation ${sample_id} contigs"}
    
    publishDir "$params.outdir/prokka", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs) 
    
    output:
    path '*.gff3', emit: gff
    path '*.gbk', emit: gbk
    
    script:
    
    """
    prokka --outdir ${sample_id}_prokka --prefix ${sample_id} ${contigs}
    """
}


process ROARY {
    label 'roary'
    tag "Roary pangenome analysis"

    publishDir "${params.outdir}/roary", mode: 'copy'

    input:
    // Accept a LIST of staged GFFs (emitted by `collect()`)
    path annotations

    output:
    path "roary_output"

    script:
    """
    roary -e --mafft -p ${task.cpus} -f roary_output ${annotations.join(' ')}
    """
}

process PHYLOGENY {
    label 'phylogeny'
    tag {"Phylogenetic tree construction"}
    
    publishDir "$params.outdir/phylogeny", mode: 'copy'

    input:
    path(gene_alignment) from ROARY.out.map { file("${it}/core_gene_alignment.aln") }
    
    output:
    path ("phylogenetic_tree.nwk")
    
    script:
    
    """
    FastTree -nt ${gene_alignment} > phylogenetic_tree.nwk
    """
}   

process KRAKEN2 {
    label 'kraken2'
    tag {"Kraken2 classification ${sample_id} reads"}
    
    publishDir "$params.outdir/kraken2", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) 
    
    output:
    path ("${sample_id}_kraken2_report.txt"), path("${sample_id}_kraken2_output.txt")
    
    script:
    
    """
    kraken2 --db ${params.kraken2_db} --threads ${task.cpus} --report ${sample_id}_kraken2_report.txt --output ${sample_id}_kraken2_output.txt ${reads}
    """
}   

process COUNT_BASES_CALLED {
    label 'count_bases'
    tag {"Count bases called ${sample_id} reads"}
    
    publishDir "$params.outdir/base_counts", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) 
    
    output:
    path ("${sample_id}_base_count.txt")
    
    script:
    
    """
    seqtk fqchk ${reads} | grep 'bases' > ${sample_id}_base_count.txt
    """
}

process GENOME_DEPTH {
    label 'genome_depth'
    tag {"Genome depth ${sample_id} reads"}
    
    publishDir "$params.outdir/genome_depth", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(contigs) 
    
    output:
    path ("${sample_id}_genome_depth.txt")
    
    script:
    
    """
    minimap2 -a -x map-ont ${contigs} ${reads} | samtools sort -o ${sample_id}_sorted.bam
    samtools index ${sample_id}_sorted.bam
    samtools depth ${sample_id}_sorted.bam > ${sample_id}_genome_depth.txt
    """
}

process SUMMARY_BLASTN {
    label 'summary_blastn'
    tag {"Summary BLASTN ${sample_id} reads"}
    
    publishDir "$params.outdir/summary_blastn", mode: 'copy'

    input:
    tuple val(sample_id), path(blastn) 
    
    output:
    path ("${sample_id}_summary_blastn.txt")
    
    script:
    
    """
    awk -v OFS='\\t' '{print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12}' ${blastn} > ${sample_id}_summary_blastn.txt
    """
}

process IQTREE {
    label 'iqtree'
    tag {"IQ-TREE phylogenetic analysis"}
    
    publishDir "$params.outdir/iqtree", mode: 'copy'

    input:
    path(gene_alignment) from ROARY.out.map { file("${it}/core_gene_alignment.aln") }
    
    output:
    path ("iqtree_output")
    
    script:
    
    """
    iqtree -s ${gene_alignment} -m GTR+G -bb 1000 -nt AUTO -pre iqtree_output
    """
}

process RAXML {
    label 'raxml'
    tag {"RAxML phylogenetic analysis"}
    
    publishDir "$params.outdir/raxml", mode: 'copy'

    input:
    path(gene_alignment) from ROARY.out.map { d -> file("${d}/core_gene_alignment.aln") }
    
    output:
    path ("raxml_output")
    
    script:
    
    """
    raxmlHPC -s ${gene_alignment} -n raxml_output -m GTRGAMMA -p 12345 -# 100 -w \$(pwd)/raxml_output
    """
}

