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
include { COUNT_BASES_CALLED } from './modules/nanopore.nf'
include { GENOME_DEPTH } from './modules/nanopore.nf'
include { ASSEMBLY_STATS; ASSEMBLY_QUALITY; GENOME_DEPTH } from './modules/nanopore.nf'
include { RAXML } from './modules/nanopore.nf'
include { IQTREE } from './modules/nanopore.nf'
include { GTDBTK } from './modules/nanopore.nf'
include { AMR_ABRFORMAT } from './modules/nanopore.nf'
include { PROKKA } from './modules/nanopore.nf'
include { ROARY } from './modules/nanopore.nf'
include { RAWFASTQC_SINGLE; CLEANFASTQC_SINGLE; MULTIQC_READS } from './modules/nanopore.nf'
include { FILTLONG } from './modules/nanopore.nf'
include { FASTP_SINGLE } from './modules/nanopore.nf'
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
           .set{fast5s}
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
           .set{refFasta}
       main:
       RAWFASTQC_SINGLE(reads)
       //FASTP_SINGLE(reads)
       FILTLONG(reads)
       CLEANFASTQC_SINGLE(FASTP_SINGLE.out.cat_fastq)
       MULTIQC_READS(CLEANFASTQC_SINGLE.out.collect())
       ASSEMBLY(FILTLONG.out.reads)
       QUAST_FROM_READS(ASSEMBLY.out.assembly)
       MULTIQC_CONTIGS(QUAST_FROM_READS.out.quast_dir.collect())
       GTDBTK(ASSEMBLY.out.assembly)
       AMR_ABRFORMAT(ASSEMBLY.out.assembly)
       PROKKA(ASSEMBLY.out.assembly)
       ROARY(PROKKA.out.gff.collect())
       RAXML(ROARY.out.core_gene_alignment)
       IQTREE(RAXML.out.raxml_bestTree)
}

