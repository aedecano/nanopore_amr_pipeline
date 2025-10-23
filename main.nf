// ================== Nextflow Pipeline for Nanopore AMR and Pangenome Analysis ==================

nextflow.enable.dsl = 2

// Load modules
include { NANOPLOT_RAW; FILTLONG; NANOPLOT_FILT; FLYE; QUAST; BANDAGE_IMAGE; ABRICATE; ABRICATE_TAG as ABRICATE_TAG_AMR; ABRICATE_TAG as ABRICATE_TAG_PLM; MERGE_ABRICATE as MERGE_ABRICATE_AMR; MERGE_ABRICATE as MERGE_ABRICATE_PLM; COMBINE_ABRICATE_RESULTS; PLOT_SUMMARIZE_ABRICATE; MLST_CONTIGS_PS; MERGE_MLST; BAKTA; MERGE_GFF_FASTA; PANAROO; SELECT_PANAROO_GML; PLOT_PANAROO_GML; RAXML_NG; IQTREE2; MULTIQC; KRAKEN2; GTDBTK_CLASSIFY } from './modules/nanopore.nf'


// ---- Params defaults (DSL2-safe) ----
params.help       = params.help ?: false
params.reads      = params.reads ?: null
params.assemblies = params.assemblies ?: null
params.outdir     = params.outdir ?: 'results'

// Optional help guard (leave out if you don’t print HELP elsewhere)
// if (params.help) { log.info "Use --reads or --assemblies with an -entry"; System.exit(0) }

def HELP = """
USAGE
  # Reads workflow
  nextflow run ./main.nf -entry assembly_amr_pangenome \\
    --reads '../ONT_Kpn_FQ/*.fastq.gz' --outdir out

  # Assemblies workflows
  nextflow run ./main.nf -entry amr_annotation_from_assemblies \\
    --assemblies 'path/*/*.contigs.fasta' --outdir out

  nextflow run ./main.nf -entry amr_pangenome_from_assemblies \\
    --assemblies 'path/*/*.contigs.fasta' --outdir out
"""

if (params.help) { log.info HELP; System.exit(0) }


// ---------------- Workflows ----------------

// -------- Workflow: assembly_amr_pangenome --------
workflow assembly_amr_pangenome {

  if (!params.reads) error "Provide --reads '<glob>' for this workflow"
log.info "reads pattern: ${params.reads}"

reads = Channel
  .fromPath(params.reads, checkIfExists: true)
  .ifEmpty { error "No files matched: ${params.reads}" }
  .map { f ->
    def sid = f.name
      .replaceFirst(/\.fastq\.gz$/, '')
      .replaceFirst(/\.fq\.gz$/,    '')
      .replaceFirst(/\.fastq$/,     '')
      .replaceFirst(/\.fq$/,        '')
    tuple(sid, f)
  }

  main:
    // QC -> Filter -> QC
    np_raw   = NANOPLOT_RAW(reads)
    filt     = FILTLONG(reads)
    np_filt  = NANOPLOT_FILT(filt.reads)

    // Assemble
    flye     = FLYE(filt.reads)

    // Assembly evaluation
    quast = QUAST(flye.assembly)

    // Bandage exports (runs only for samples where a GFA exists)
    bandage = BANDAGE_IMAGE(flye.graph)

    // ABRicate (per sample), tag, merge
    ab              = ABRICATE(flye.assembly)
    tag_amr         = ABRICATE_TAG_AMR(ab.amr_tsv, 'amr')
    tag_plm         = ABRICATE_TAG_PLM(ab.plm_tsv, 'plm')
    merged_amr      = MERGE_ABRICATE_AMR(tag_amr.tagged.collect(), 'amr')
    merged_plm      = MERGE_ABRICATE_PLM(tag_plm.tagged.collect(), 'plm')
    abri_all        = ab.amr_tsv.mix(ab.plm_tsv).collect()
    abri_merged     = merged_amr.merged.mix(merged_plm.merged).collect()
    combined        = COMBINE_ABRICATE_RESULTS(merged_amr.merged, merged_plm.merged)
    abri_summ_plots = PLOT_SUMMARIZE_ABRICATE(combined.combined)

    // MultiQC inputs
    mqc_inputs = np_raw.report_dir.map { _, d -> d }
                .mix( quast.report_dir.map { _, d -> d } )
                .mix( np_filt.report_dir.map { _, d -> d } )
                .mix( tag_amr.tagged )
                .mix( tag_plm.tagged )
                .collect()
    mqc = MULTIQC(mqc_inputs)

    // MLST (per sample, merge then summarise)
    mlst       = MLST_CONTIGS_PS(flye.assembly)
    merge_mlst = MERGE_MLST( mlst.mlst_tsv.map{ it[1] }.collect() )

    // Bakta (per sample)
    bakta    = BAKTA(flye.assembly)

    // === CORRECTED PANAROO INTEGRATION ===
    // Join Bakta GFF with assembly on sample_id
    paired_annot = bakta.gff.join(flye.assembly)
    
    // Merge GFF + FASTA for each sample
    merged_gffs = MERGE_GFF_FASTA(paired_annot)
    
    // Collect all merged GFFs and run Panaroo
    gff_list = merged_gffs.merged.map { sid, gff -> gff }.collect()
    pana = PANAROO(gff_list)

    // Plot Panaroo results
    pana_gml = SELECT_PANAROO_GML(pana.dir)
    pana_plots = PLOT_PANAROO_GML(pana_gml)

    // Optional: Trees from core alignment
    // raxml = RAXML_NG(pana.fasta)
    // iq = IQTREE2(pana.fasta)

    
  emit:
    nanoplot_raw_dirs    = np_raw.report_dir
    nanoplot_filt_dirs   = np_filt.report_dir
    filtered_reads       = filt.reads
    assemblies           = flye.assembly
    gfa_files            = flye.graph
    asm_info_files       = flye.info
    quast_dir            = quast.report_dir
    quast_tsv            = quast.report_tsv
    quast_txt            = quast.report_txt
    bandage_png          = bandage.png
    bandage_svg          = bandage.svg
    bandage_nodes_table  = bandage.info
    abricate_amr_per_sample  = ab.amr_tsv
    abricate_amr_summary     = ab.amr_summary
    abricate_plm_per_sample  = ab.plm_tsv
    abricate_plm_summary     = ab.plm_summary
    abricate_all_per_sample  = abri_all
    abricate_merged_tsv      = abri_merged
    abricate_combined        = combined.combined
    abricate_summary_plots   = abri_summ_plots.abricate_summary_plots
    per_sample_summary       = abri_summ_plots.per_sample_summary
    mlst_tsv             = mlst.mlst_tsv
    mlst_merged          = merge_mlst.mlst_merged
    mlst_summary         = merge_mlst.mlst_summary  
    bakta_gff            = bakta.gff
    bakta_tsv            = bakta.tsv
    bakta_faa            = bakta.faa
    bakta_ffn            = bakta.ffn
    bakta_txt            = bakta.txt
    bakta_json           = bakta.json
    bakta_embl           = bakta.embl
    bakta_gbff           = bakta.gbff
    bakta_inference_tsv  = bakta.inference_tsv
    bakta_svg            = bakta.svg
    bakta_png            = bakta.png
    bakta_log            = bakta.log
    panaroo_dir          = pana.dir
    core_alignment       = pana.fasta
    roary_like           = pana.csv
    graph_gml            = pana.gml
    panaroo_plots        = pana_plots.all_plots
    panaroo_plot_png     = pana_plots.plot
    multiqc_report       = mqc.report
}


// -------- Workflow: taxonomy_from_reads --------
workflow taxonomy_from_reads {

  if (!params.reads) error "Provide --reads '<glob>' for this workflow"
  log.info "reads pattern: ${params.reads}"

  reads = Channel
  .fromPath(params.reads, checkIfExists: true)
  .ifEmpty { error "No files matched: ${params.reads}" }
  .map { f ->
    def sid = f.name
      .replaceFirst(/\.fastq\.gz$/, '')
      .replaceFirst(/\.fq\.gz$/,    '')
      .replaceFirst(/\.fastq$/,     '')
      .replaceFirst(/\.fq$/,        '')
    tuple(sid, f)
  }

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


// -------- Workflow: amr_annotation_from_assemblies --------
// Takes assemblies, runs QUAST, ABRicate, MLST, and Bakta (no pangenome)
workflow amr_annotation_from_assemblies {

  if (!params.assemblies) error "Provide --assemblies '<glob>' for this workflow" 
  log.info "assemblies pattern: ${params.assemblies}"

  assemblies = Channel
  .fromPath(params.assemblies, checkIfExists: true)
  .ifEmpty { error "No assemblies matched: ${params.assemblies}" }
  .map { f ->
    def sid = f.getBaseName()
      .replaceFirst(/\.fa(sta)?(\.gz)?$/, '')
      .replaceFirst(/\.fna(\.gz)?$/,      '')
    tuple(sid, f)
  }

  main:

    asm = assemblies
    // Assembly evaluation
    quast = QUAST(asm)

    // Summarise Quast results with MultiQC
    mqc_quast = MULTIQC(quast.report_dir.map { _, d -> d }.collect())
    
    // ABRicate (per sample), tag, merge
    ab          = ABRICATE(asm)
    tag_amr     = ABRICATE_TAG_AMR(ab.amr_tsv, 'amr')
    tag_plm     = ABRICATE_TAG_PLM(ab.plm_tsv, 'plm')
    merged_amr  = MERGE_ABRICATE_AMR(tag_amr.tagged.collect(), 'amr')
    merged_plm  = MERGE_ABRICATE_PLM(tag_plm.tagged.collect(), 'plm')
    abri_all    = ab.amr_tsv.mix(ab.plm_tsv).collect()
    abri_merged = merged_amr.merged.mix(merged_plm.merged).collect()
    combined    = COMBINE_ABRICATE_RESULTS(merged_amr.merged, merged_plm.merged)
    abri_summ_plots = PLOT_SUMMARIZE_ABRICATE(combined.combined)

    // MLST (per sample, merge then summarise)
    mlst       = MLST_CONTIGS_PS(asm)
    merge_mlst = MERGE_MLST( mlst.mlst_tsv.map{ it[1] }.collect() )

    // Bakta (per sample) - FINAL STEP
    //bakta = BAKTA(asm)

  emit:
    quast_dir                = quast.report_dir
    quast_tsv                = quast.report_tsv
    quast_txt                = quast.report_txt
    multiqc_report           = mqc_quast.report
    abricate_amr_per_sample  = ab.amr_tsv
    abricate_amr_summary     = ab.amr_summary
    abricate_plm_per_sample  = ab.plm_tsv
    abricate_plm_summary     = ab.plm_summary
    abricate_all_per_sample  = abri_all
    abricate_merged_tsv      = abri_merged
    abricate_combined        = combined.combined
    abricate_summary_plots   = abri_summ_plots.abricate_summary_plots
    per_sample_summary       = abri_summ_plots.per_sample_summary
    mlst_tsv                 = mlst.mlst_tsv
    mlst_merged              = merge_mlst.mlst_merged
    mlst_summary             = merge_mlst.mlst_summary  
    //bakta_gff                = bakta.gff
    //bakta_tsv                = bakta.tsv
    //bakta_faa                = bakta.faa
    //bakta_ffn                = bakta.ffn
    //bakta_txt                = bakta.txt
    //bakta_json               = bakta.json
    //bakta_embl               = bakta.embl
    //bakta_gbff               = bakta.gbff
    //bakta_inference_tsv      = bakta.inference_tsv
    //bakta_svg                = bakta.svg
    //bakta_png                = bakta.png
    //bakta_log                = bakta.log
}


// -------- Workflow: amr_pangenome_from_assemblies --------
// Takes assemblies, runs QUAST, ABRicate, MLST, Bakta, and Panaroo pangenome
workflow amr_pangenome_from_assemblies {
   
  if (!params.assemblies) error "Provide --assemblies '<glob>' for this workflow"
  log.info "assemblies pattern: ${params.assemblies}"

  asm = Channel
  .fromPath(params.assemblies, checkIfExists: true)
  .ifEmpty { error "No assemblies matched: ${params.assemblies}" }
  .map { f ->
    def sid = f.getBaseName()
      .replaceFirst(/\.fa(sta)?(\.gz)?$/, '')
      .replaceFirst(/\.fna(\.gz)?$/,      '')
    tuple(sid, f)
  }

  main:
    
    // Assembly evaluation
    //quast = QUAST(asm)

    // Summarise Quast results with MultiQC
    //mqc_quast = MULTIQC(quast.report_dir.map { _, d -> d }.collect())
    
    // ABRicate (per sample), tag, merge
    //ab          = ABRICATE(asm)
    //tag_amr     = ABRICATE_TAG_AMR(ab.amr_tsv, 'amr')
    //tag_plm     = ABRICATE_TAG_PLM(ab.plm_tsv, 'plm')
    //merged_amr  = MERGE_ABRICATE_AMR(tag_amr.tagged.collect(), 'amr')
    //merged_plm  = MERGE_ABRICATE_PLM(tag_plm.tagged.collect(), 'plm')
    //abri_all    = ab.amr_tsv.mix(ab.plm_tsv).collect()
    //abri_merged = merged_amr.merged.mix(merged_plm.merged).collect()
    //combined    = COMBINE_ABRICATE_RESULTS(merged_amr.merged, merged_plm.merged)
    //abri_summ_plots = PLOT_SUMMARIZE_ABRICATE(combined.combined)

    // MLST (per sample, merge then summarise)
    //mlst       = MLST_CONTIGS_PS(asm)
    //merge_mlst = MERGE_MLST( mlst.mlst_tsv.map{ it[1] }.collect() )

    // Prepare assemblies channel keyed by sample_id for Panaroo input downstream
    assemblies_keyed = asm


    // Bakta (per sample)
    bakta = BAKTA(assemblies_keyed)


    // Join Bakta GFF with assembly on sample_id
    annot_with_asm = bakta.gff.join(assemblies_keyed)
    
    // Merge GFF + FASTA for each sample and save the list of merged GFFs
    merged = MERGE_GFF_FASTA(annot_with_asm)
    merged_gff_list = merged
      .map { sid, merged_gff -> merged_gff }
      .collectFile(name: 'merged_gff_list.txt')
    
    // Run Panaroo on the merged list of GFFs
    pana = PANAROO(merged_gff_list)

    // Plot Panaroo results
    //pana_gml = SELECT_PANAROO_GML(pana.dir)
    //pana_plots = PLOT_PANAROO_GML(pana_gml)

    // Optional: Trees from core alignment
    // raxml = RAXML_NG(pana.fasta)
    // iq = IQTREE2(pana.fasta)

  emit:
    //quast_dir                = quast.report_dir
    //quast_tsv                = quast.report_tsv
    //quast_txt                = quast.report_txt
    //multiqc_report           = mqc_quast.report
    //abricate_amr_per_sample  = ab.amr_tsv
    //abricate_amr_summary     = ab.amr_summary
    //abricate_plm_per_sample  = ab.plm_tsv
    //abricate_plm_summary     = ab.plm_summary
    //abricate_all_per_sample  = abri_all
    //abricate_merged_tsv      = abri_merged
    //abricate_combined        = combined.combined
    //abricate_summary_plots   = abri_summ_plots.abricate_summary_plots
    //per_sample_summary       = abri_summ_plots.per_sample_summary
   //mlst_tsv                 = mlst.mlst_tsv
    //mlst_merged              = merge_mlst.mlst_merged
    //mlst_summary             = merge_mlst.mlst_summary  
    emit:
    bakta_gff                = bakta.gff
    bakta_tsv                = bakta.tsv
    bakta_faa                = bakta.faa
    bakta_ffn                = bakta.ffn
    bakta_txt                = bakta.txt
    bakta_json               = bakta.json
    bakta_embl               = bakta.embl
    bakta_gbff               = bakta.gbff
    bakta_inference_tsv      = bakta.inference_tsv
    bakta_svg                = bakta.svg
    bakta_png                = bakta.png
    bakta_log                = bakta.log
    panaroo_dir              = pana.dir
    core_aln                 = pana.core_aln
    roary_like               = pana.roary_like
    graph_gml                = pana.graph_gml     
    graph_gml_alt            = pana.graph_gml_alt
    //panaroo_plots            = pana_plots.all_plots
    //panaroo_plot_png         = pana_plots.plot
}

// -------- Workflow: assembly_amr_annotation --------

// Takes nanopore reads, does QC, filtering, assembly, evaluation,
// AMR/plasmid detection, MLST, and Bakta annotation (no pangenome analysis)

workflow assembly_amr_annotation {
  if (!params.reads) error "Provide --reads '<glob>' for this workflow"
  log.info "reads pattern: ${params.reads}"

  reads = Channel
  .fromPath(params.reads, checkIfExists: true)
  .ifEmpty { error "No files matched: ${params.reads}" }
  .map { f ->
    def sid = f.name
      .replaceFirst(/\.fastq\.gz$/, '')
      .replaceFirst(/\.fq\.gz$/,    '')
      .replaceFirst(/\.fastq$/,     '')
      .replaceFirst(/\.fq$/,        '')
    tuple(sid, f)
  }

  main:
    // QC -> Filter -> QC
    np_raw   = NANOPLOT_RAW(reads)
    filt     = FILTLONG(reads)
    np_filt  = NANOPLOT_FILT(filt.reads)

    // Assemble
    flye     = FLYE(filt.reads)

    // Assembly evaluation
    quast = QUAST(flye.assembly)

    // Bandage exports (runs only for samples where a GFA exists)
    bandage = BANDAGE_IMAGE(flye.graph)

    // ABRicate (per sample), tag, merge
    ab              = ABRICATE(flye.assembly)
    tag_amr         = ABRICATE_TAG_AMR(ab.amr_tsv, 'amr')
    tag_plm         = ABRICATE_TAG_PLM(ab.plm_tsv, 'plm')
    merged_amr      = MERGE_ABRICATE_AMR(tag_amr.tagged.collect(), 'amr')
    merged_plm      = MERGE_ABRICATE_PLM(tag_plm.tagged.collect(), 'plm')
    abri_all        = ab.amr_tsv.mix(ab.plm_tsv).collect()
    abri_merged     = merged_amr.merged.mix(merged_plm.merged).collect()
    combined        = COMBINE_ABRICATE_RESULTS(merged_amr.merged, merged_plm.merged)
    abri_summ_plots = PLOT_SUMMARIZE_ABRICATE(combined.combined)

    // MultiQC inputs
    mqc_inputs = np_raw.report_dir.map { _, d -> d }
                .mix( quast.report_dir.map { _, d -> d } )
                .mix( np_filt.report_dir.map { _, d -> d } )
                .mix( tag_amr.tagged )
                .mix( tag_plm.tagged )
                .collect()
    mqc = MULTIQC(mqc_inputs)

    // MLST (per sample, merge then summarise)
    mlst       = MLST_CONTIGS_PS(flye.assembly)
    merge_mlst = MERGE_MLST( mlst.mlst_tsv.map{ it[1] }.collect() )

    // Bakta (per sample) - FINAL STEP
    bakta    = BAKTA(flye.assembly)

  emit:
    nanoplot_raw_dirs        = np_raw.report_dir
    nanoplot_filt_dirs       = np_filt.report_dir
    filtered_reads           = filt.reads
    assemblies               = flye.assembly
    gfa_files                = flye.graph
    asm_info_files           = flye.info
    quast_dir                = quast.report_dir
    quast_tsv                = quast.report_tsv
    quast_txt                = quast.report_txt
    bandage_png              = bandage.png
    bandage_svg              = bandage.svg
    bandage_nodes_table      = bandage.info
    abricate_amr_per_sample  = ab.amr_tsv
    abricate_amr_summary     = ab.amr_summary
    abricate_plm_per_sample  = ab.plm_tsv
    abricate_plm_summary     = ab.plm_summary
    abricate_all_per_sample  = abri_all
    abricate_merged_tsv      = abri_merged
    abricate_combined        = combined.combined
    abricate_summary_plots   = abri_summ_plots.abricate_summary_plots
    per_sample_summary       = abri_summ_plots.per_sample_summary
    mlst_tsv                 = mlst.mlst_tsv
    mlst_merged              = merge_mlst.mlst_merged
    mlst_summary             = merge_mlst.mlst_summary  
    bakta_gff                = bakta.gff
    bakta_tsv                = bakta.tsv
    bakta_faa                = bakta.faa
    bakta_ffn                = bakta.ffn
    bakta_txt                = bakta.txt
    bakta_json               = bakta.json
    bakta_embl               = bakta.embl
    bakta_gbff               = bakta.gbff
    bakta_inference_tsv      = bakta.inference_tsv
    bakta_svg                = bakta.svg
    bakta_png                = bakta.png
    bakta_log                = bakta.log
    multiqc_report           = mqc.report
}


// Default workflow entry point
workflow {
  assembly_amr_pangenome()
}