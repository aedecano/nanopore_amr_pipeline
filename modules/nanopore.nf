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
    tag {"Dorado basecalling ${uuid} reads"}
    
    publishDir "$params.outdir/basecalled", mode: 'copy'

    input:
    tuple val(uuid), path(reads) 
    
    output:
    tuple val(uuid), path("${uuid}_fastq/${uuid}_1.fastq.gz"), path("${uuid}_fastq/${uuid}_2.fastq.gz")
    
    script:
    
    """
    dorado dna_r9.4.1 -i ${reads} -o ${uuid}_fastq -t ${task.cpus} --device auto --bam_out ${uuid}_fastq/${uuid}.bam
    """
}

process BASECALL_GUPPY {
    label 'guppy basecall'
    tag {"Guppy basecalling ${uuid} reads"}
    
    publishDir "$params.outdir/basecalled", mode: 'copy'

    input:
    tuple val(uuid), path(reads) 
    
    output:
    tuple val(uuid), path("${uuid}_fastq/${uuid}_1.fastq.gz"), path("${uuid}_fastq/${uuid}_2.fastq.gz")
    
    script:
    
    """
    guppy_basecaller -i ${reads} -s ${uuid}_fastq --flowcell FLO-MIN106 --kit SQK-LSK109 --num_callers ${task.cpus} --cpu_threads_per_caller 4 --barcode_kits EXP-NBD104 --trim_barcodes
    """
}


process RAWFASTQC_SINGLE {
    label 'fastqc'
    tag {"FastQC raw ${uuid} reads"}
    
    publishDir "$params.outdir/raw_fastqc", mode: 'copy'

    input:
    tuple val(uuid), path(reads) 
    
    output:
    path ("${uuid}_fastqc/${uuid}.txt")
    
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
    tag {"Filtlong ${uuid} reads"}
    
    publishDir "$params.outdir/cleaned_reads", mode: 'copy'

    input:
    tuple val(uuid), path(reads) 
    
    output:
    tuple val(uuid), path("${uuid}_filtlong.fastq")
    
    script:
    
    """
    filtlong --min_length 1000 --keep_percent 90 ${reads} > ${uuid}_filtlong.fastq
    """
}

process CLEANFASTQC_SINGLE {
    label 'fastqc'
    tag {"FastQC cleaned ${uuid} reads"}
    
    publishDir "$params.outdir/cleaned_fastqc", mode: 'copy'

    input:
    tuple val(uuid), path(reads) 
    
    output:
    path ("${uuid}_clean_fastqc/${uuid}.txt")
    
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

process ASSEMBLY {
    label 'assembly'
    tag {"Assembly ${uuid} reads"}
    
    publishDir "$params.outdir/assemblies", mode: 'copy'

    input:
    tuple val(uuid), path(reads) 
    
    output:
    tuple val(uuid), path("${uuid}_assembly/${uuid}_contigs.fasta")
    
    script:
    
    """
    flye --nano-raw ${reads} --out-dir ${uuid}_assembly --threads ${task.cpus} --genome-size 5m
    """
}

process QUAST_FROM_READS {
    label 'quast'
    tag {"QUAST from ${uuid} reads"}
    
    publishDir "$params.outdir/quast_reads", mode: 'copy'

    input:
    tuple val(uuid), path(reads) 
    
    output:
    path ("${uuid}_quast_report")
    
    script:
    
    """
    quast.py ${reads} -o ${uuid}_quast_report -t ${task.cpus}
    """
}   

process QUAST_FROM_CONTIGS {
    label 'quast'
    tag {"QUAST from ${uuid} contigs"}
    
    publishDir "$params.outdir/quast_contigs", mode: 'copy'

    input:
    tuple val(uuid), path(contigs) 
    
    output:
    path ("${uuid}_quast_report")
    
    script:
    
    """
    quast.py ${contigs} -o ${uuid}_quast_report -t ${task.cpus}
    """
}   

process MEDAKA {
    label 'medaka'
    tag {"Medaka polishing ${uuid} reads"}
    
    publishDir "$params.outdir/polished_contigs", mode: 'copy'

    input:
    tuple val(uuid), path(reads), path(contigs) 
    
    output:
    tuple val(uuid), path("${uuid}_polished/${uuid}_polished.fasta")
    
    script:
    
    """
    medaka_consensus -i ${reads} -d ${contigs} -o ${uuid}_polished -t ${task.cpus} -m r941_min_high_g360
    """
}

process NANOPLOTYPER {
    label 'nanoplot'
    tag {"NanoPlot ${uuid} reads"}
    
    publishDir "$params.outdir/nanoplot", mode: 'copy'

    input:
    tuple val(uuid), path(reads) 
    
    output:
    path ("${uuid}_nanoplot")
    
    script:
    
    """
    NanoPlot --fastq ${reads} -o ${uuid}_nanoplot --threads ${task.cpus}
    """
}

process GTDBTYPER {
    label 'gtdbtyper'
    tag {"GTDB-Tk classification ${uuid} contigs"}
    
    publishDir "$params.outdir/gtdbtk", mode: 'copy'

    input:
    tuple val(uuid), path(contigs) 
    
    output:
    path ("${uuid}_gtdbtk_report.txt")
    
    script:
    
    """
    gtdbtk classify_wf --genome_dir . --out_dir ${uuid}_gtdbtk --cpus ${task.cpus}
    cp ${uuid}_gtdbtk/gtdbtk.bac120.summary.tsv ${uuid}_gtdbtk_report.txt
    """
}


process AMR_PLM_FROM_CONTIGS {
    label 'amr_plm'
    tag {"AMR Plasmid from ${uuid} contigs"}
    
    publishDir "$params.outdir/amr_plasmids_contigs", mode: 'copy'

    input:
    tuple val(uuid), path(contigs) 
    
    output:
    path ("${uuid}_plasmid_report.txt")
    
    script:
    
    """
    plasmidfinder.py -i ${contigs} -o ${uuid}_plasmid_report -p /path/to/plasmidfinder_db -t 0.95 -l 0.6
    """
}

process AMR_ABRFORMAT {
    label 'amr_abrformat'
    tag {"AMR ABRicate from ${uuid} contigs"}
    
    publishDir "$params.outdir/amr_abricate_contigs", mode: 'copy'

    input:
    tuple val(uuid), path(contigs) 
    
    output:
    path ("${uuid}_abricate_report.txt")
    
    script:
    
    """
    abricate --db ncbi ${contigs} > ${uuid}_abricate_report.txt
    """
}

process MOBTYPER {
    label 'mobtyper'
    tag {"MOB-typer from ${uuid} contigs"}
    
    publishDir "$params.outdir/mobtyper_contigs", mode: 'copy'

    input:
    tuple val(uuid), path(contigs) 
    
    output:
    path ("${uuid}_mobtyper_report.txt")
    
    script:
    
    """
    mob_typer -i ${contigs} -o ${uuid}_mobtyper_report
    """
}

process AMRFINDERPLUS {
    label 'amrfinderplus'
    tag {"AMR Finder Plus from ${uuid} contigs"}
    
    publishDir "$params.outdir/amrfinderplus_contigs", mode: 'copy'

    input:
    tuple val(uuid), path(contigs) 
    
    output:
    path ("${uuid}_amrfinderplus_report.txt")
    
    script:
    
    """
    amrfinder -n ${contigs} -o ${uuid}_amrfinderplus_report.txt --organism Clostridioides_difficile --plus
    """
}

process MLST_FROM_CONTIGS {
    label 'mlst'
    tag {"MLST from ${uuid} contigs"}
    
    publishDir "$params.outdir/mlst_contigs", mode: 'copy'

    input:
    tuple val(uuid), path(contigs) 
    
    output:
    path ("${uuid}_mlst_report.txt")
    
    script:
    
    """
    mlst ${contigs} > ${uuid}_mlst_report.txt
    """
}

process PROKKA {
    label 'prokka'
    tag {"Prokka annotation ${uuid} contigs"}
    
    publishDir "$params.outdir/prokka_annotations", mode: 'copy'

    input:
    tuple val(uuid), path(contigs) 
    
    output:
    path ("${uuid}_prokka")
    
    script:
    
    """
    prokka --outdir ${uuid}_prokka --prefix ${uuid} --cpus ${task.cpus} ${contigs}
    """
}

process ROARY {
    label 'roary'
    tag {"Roary pangenome analysis"}
    
    publishDir "$params.outdir/roary", mode: 'copy'

    input:
    path(annotations) from PROKKA.out.collect()
    
    output:
    path ("roary_output")
    
    script:
    
    """
    roary -e --mafft -p ${task.cpus} -f roary_output ${annotations}
    """
}

process PHYLOGENY {
    label 'phylogeny'
    tag {"Phylogenetic tree construction"}
    
    publishDir "$params.outdir/phylogeny", mode: 'copy'

    input:
    path(gene_alignment) from ROARY.out.map { it -> file("${it}/core_gene_alignment.aln") }
    
    output:
    path ("phylogenetic_tree.nwk")
    
    script:
    
    """
    FastTree -nt ${gene_alignment} > phylogenetic_tree.nwk
    """
}   

process KRAKEN2 {
    label 'kraken2'
    tag {"Kraken2 classification ${uuid} reads"}
    
    publishDir "$params.outdir/kraken2", mode: 'copy'

    input:
    tuple val(uuid), path(reads) 
    
    output:
    path ("${uuid}_kraken2_report.txt"), path("${uuid}_kraken2_output.txt")
    
    script:
    
    """
    kraken2 --db ${params.kraken2_db} --threads ${task.cpus} --report ${uuid}_kraken2_report.txt --output ${uuid}_kraken2_output.txt ${reads}
    """
}   

process COUNT_BASES_CALLED {
    label 'count_bases'
    tag {"Count bases called ${uuid} reads"}
    
    publishDir "$params.outdir/base_counts", mode: 'copy'

    input:
    tuple val(uuid), path(reads) 
    
    output:
    path ("${uuid}_base_count.txt")
    
    script:
    
    """
    seqtk fqchk ${reads} | grep 'bases' > ${uuid}_base_count.txt
    """
}

process GENOME_DEPTH {
    label 'genome_depth'
    tag {"Genome depth ${uuid} reads"}
    
    publishDir "$params.outdir/genome_depth", mode: 'copy'

    input:
    tuple val(uuid), path(reads), path(contigs) 
    
    output:
    path ("${uuid}_genome_depth.txt")
    
    script:
    
    """
    minimap2 -a -x map-ont ${contigs} ${reads} | samtools sort -o ${uuid}_sorted.bam
    samtools index ${uuid}_sorted.bam
    samtools depth ${uuid}_sorted.bam > ${uuid}_genome_depth.txt
    """
}

process SUMMARY_BLASTN {
    label 'summary_blastn'
    tag {"Summary BLASTN ${uuid} reads"}
    
    publishDir "$params.outdir/summary_blastn", mode: 'copy'

    input:
    tuple val(uuid), path(blastn) 
    
    output:
    path ("${uuid}_summary_blastn.txt")
    
    script:
    
    """
    awk -v OFS='\\t' '{print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12}' ${blastn} > ${uuid}_summary_blastn.txt
    """
}

process IQTREE {
    label 'iqtree'
    tag {"IQ-TREE phylogenetic analysis"}
    
    publishDir "$params.outdir/iqtree", mode: 'copy'

    input:
    path(gene_alignment) from ROARY.out.map { it -> file("${it}/core_gene_alignment.aln") }
    
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
    path(gene_alignment) from ROARY.out.map { it -> file("${it}/core_gene_alignment.aln") }
    
    output:
    path ("raxml_output")
    
    script:
    
    """
    raxmlHPC -s ${gene_alignment} -n raxml_output -m GTRGAMMA -p 12345 -# 100 -w $(pwd)/raxml_output
    """
}

