nextflow.enable.dsl = 2

include { NANOPLOT_RAW; FILTLONG; NANOPLOT_FILT; FLYE; QUAST; BANDAGE_IMAGE; ABRICATE; ABRICATE_TAG; MERGE_ABRICATE; BAKTA; PANAROO; RAXML_NG; IQTREE2; MULTIQC; KRAKEN2; GTDBTK_CLASSIFY } from './modules/nanopore.nf'

// -------- Params --------
params.kraken2_db = params.kraken2_db ?: null
params.gtdbtk_db  = params.gtdbtk_db  ?: null
params.reads  = params.reads ?: null
params.outdir = params.outdir ?: 'results'

// -------- Channels --------
Channel
  .fromPath(params.reads, checkIfExists: true)
  .ifEmpty { error "No reads found for: ${params.reads}" }
  .map { f ->
    def sid = f.name
      .replaceFirst(/\.fastq\.gz$/, '')
      .replaceFirst(/\.fq\.gz$/, '')
      .replaceFirst(/\.fastq$/, '')
      .replaceFirst(/\.fq$/, '')
    tuple(sid, f)
  }
  .set { reads }   // (sample_id, fastq.gz)

// Optional: drop empty files early
// reads = reads.filter { sid, fq -> fq.size() > 0 }

// ---------------- Workflows ----------------
workflow assembly_amr_pangenome {

  main:
    // QC -> Filter -> QC
    np_raw   = NANOPLOT_RAW(reads)
    filt     = FILTLONG(reads)            // -> (sample_id, .filt.fastq.gz)
    np_filt  = NANOPLOT_FILT(filt.reads)

    //Taxonomy from filtered reads (contamination check)
    kraken = KRAKEN2(filt.reads)

    // Assemble
    flye     = FLYE(filt.reads)           // -> (sample_id, .contigs.fasta)

    // Assembly evaluation
    quast = QUAST(flye.assembly)

    // Bandage exports (runs only for samples where a GFA exists)
    bandage = BANDAGE_IMAGE(flye.graph)

    // ABRicate (per sample), tag, merge
    abri     = ABRICATE(flye.assembly)
    abri_tag = ABRICATE_TAG(abri.tsv)
    abri_all = abri_tag.tagged.collect()
    abri_merged = MERGE_ABRICATE(abri_all)

    // Bakta (per sample)
    //bakta    = BAKTA(flye.assembly)

    // Panaroo (needs a list of GFFs)
    //gff_list = bakta.gff.map { sid, g -> g }.collect()
    //pana     = PANAROO(gff_list)

    // Trees from core alignment
    //raxml    = RAXML_NG(pana.core_alignment)
    //iq       = IQTREE2(pana.core_alignment)

    // MultiQC inputs: NanoPlot dirs and ABRicate TSVs (MultiQC parses both)
    mqc_inputs = np_raw.report_dir.map { _, d -> d }
                .mix( quast.report_dir.map { _, d -> d } )
                .mix( np_filt.report_dir.map { _, d -> d } )
                .mix( abri.tsv.map { _, f -> f } )
                .collect()
    mqc = MULTIQC(mqc_inputs)
    
  emit:
    nanoplot_raw_dirs    = np_raw.report_dir
    nanoplot_filt_dirs   = np_filt.report_dir
    filtered_reads       = filt.reads
    kraken_report        = kraken.report
    kraken_calls         = kraken.classified
    assemblies           = flye.assembly
    gfa_files            = flye.graph
    asm_info_files       = flye.info
    quast_dir            = quast.report_dir
    quast_tsv            = quast.report_tsv
    quast_txt            = quast.report_txt
    bandage_png         = bandage.png
    bandage_svg         = bandage.svg
    bandage_nodes_table = bandage.info
    abricate_per_sample  = abri.tsv
    abricate_summary     = abri.summary
    abricate_merged_tsv  = abri_merged.merged
    //bakta_gff            = bakta.gff
    //core_alignment       = pana.core_alignment
    //raxml_tree           = raxml.besttree
    //iqtree_tree          = iq.treefile
    multiqc_report       = mqc.report
}

// default entrypoint


// -------- Workflow: taxonomy_from_reads --------
// Takes nanopore reads (tuple: sample_id, reads), runs Kraken2 on reads,
// assembles with Flye, then runs GTDB-Tk on the assembly.
workflow taxonomy_from_reads {

  take:
    reads

  main:
    // Run Kraken2 directly on raw reads
    kraken = KRAKEN2(reads)

    // Assemble reads with Flye, then classify assembly with GTDB-Tk
    flye = FLYE(reads)
    gtdb = GTDBTK_CLASSIFY(flye.assembly)

  emit:
  kraken_report = kraken.report
  kraken_calls  = kraken.classified
  assembly      = flye.assembly
  gtdb_summary  = gtdb.summary
  gtdb_class    = gtdb.classification
}

// make this the default entry with:
// workflow { taxonomy_from_reads(reads) }

// default entrypoint
workflow { assembly_amr_pangenome() }