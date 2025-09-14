#!/usr/bin/env nextflow

/* 
A pipeline for de novo assembly, specied ID and AMR annotation of Nanopore long read libraries

QC


Mapping and Variant calling
 
 
Assembly
*/

/*
#==============================================
Enable DSL2
#==============================================
*/

nextflow.enable.dsl = 2

/*
#==============================================
Modules
#==============================================
*/

include { BASECALL_DORADO; COUNT_BASES_CALLED } from './modules/nanopore.nf'
include { ASSEMBLY } from './modules/nanopore.nf'
include { KRAKEN2 } from './modules/nanopore.nf'
include { GENOME_DEPTH } from './modules/nanopore.nf'
//include { ASSEMBLY_STATS; ASSEMBLY_QUALITY; GENOME_DEPTH } from './modules/nanopore.nf'
include { RAXML } from './modules/nanopore.nf'
include { IQTREE } from './modules/nanopore.nf'
include { GTDBTYPER } from './modules/nanopore.nf'
include { AMR_ABRFORMAT } from './modules/nanopore.nf'
include { PROKKA } from './modules/nanopore.nf'
include { ROARY } from './modules/nanopore.nf'
include { RAWFASTQC_SINGLE; CLEANFASTQC_SINGLE; MULTIQC_READS } from './modules/nanopore.nf'
include { FILTLONG } from './modules/nanopore.nf'
//include { FASTP_SINGLE } from './modules/nanopore.nf'
include { QUAST_FROM_READS; MULTIQC_CONTIGS } from './modules/nanopore.nf'

/*
#==============================================
Parameters
#==============================================
*/

params.fast5 = ""
params.ref = ""
params.reads = ""
params.outdir = ""
params.contigs = ""
params.mlstdb = ""
params.prefix = ""
params.blastn = ""
params.mlst_loci = ""
params.kraken2_db= ""

workflow basecall_extractFQ {
       Channel.fromPath(params.fast5, checkIfExists: true)
           .map{it}
           //.view()
           .set{fast5}
       main:
       BASECALL_DORADO(fast5)
       BAM_TO_FASTQ(BASECALL_DORADO.out.bam)
       COUNT_BASES_CALLED(BAM_TO_FASTQ.out.fastq)
}


workflow assembly_amr_pangenome {
       Channel.fromPath(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
       Channel.fromPath(params.ref, checkIfExists:true)
           //.view()       
           .first()
           .set{ref}
       main:
       //RAWFASTQC_SINGLE(reads)
       FILTLONG(reads)
       //CLEANFASTQC_SINGLE(FILTLONG.out.reads)
       //MULTIQC_READS(RAWFASTQC_SINGLE.out.fastqc.mix(CLEANFASTQC_SINGLE.out.fastqc ))
       ASSEMBLY(FILTLONG.out.reads)
       //QUAST_FROM_READS(ASSEMBLY.out.assembly)
       //MULTIQC_CONTIGS(QUAST_FROM_READS.out.quast_dir.collect())
       //GTDBTK(ASSEMBLY.out.assembly)
       AMR_ABRFORMAT(ASSEMBLY.out.assembly)
       //PROKKA(ASSEMBLY.out.assembly)
       //ROARY(PROKKA.out.gff.collect())
       //RAXML(ROARY.out.core_gene_alignment)
       //IQTREE(RAXML.out.raxml_bestTree)
}


workflow amr_annotation {

    // Accept either a glob (e.g. ".../*.fasta") or a directory (".../Downloads")
    def pattern = params.contigs
    if (file(pattern).isDirectory()) {
        // include common FASTA extensions; adjust as you like
        pattern = "${pattern}/*.{fa,fna,fasta,fas}"
    }

    Channel
        .fromPath(pattern, checkIfExists: true)
        .map { f -> tuple(f.baseName, f) }   // (sample_id, path)
        .set { contigs }

    main:
    AMR_ABRFORMAT(contigs)
    // write a combined TSV report
    AMR_ABRFORMAT.out.abricate_report .collectFile(name: 'abricate_reports.tsv', storeDir: 'results/amr', keepHeader: true)
}


workflow kraken2_classification {
       Channel.fromPath(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
       main:
       KRAKEN2(reads, params.kraken2_db)
}

