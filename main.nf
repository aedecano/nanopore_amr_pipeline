nextflow.enable.dsl = 2

include { NANOPLOT_RAW; FILTLONG; NANOPLOT_FILT; FLYE; QUAST; BANDAGE_IMAGE; ABRICATE; ABRICATE_TAG as ABRICATE_TAG_AMR; ABRICATE_TAG as ABRICATE_TAG_PLM; MERGE_ABRICATE as MERGE_ABRICATE_AMR; MERGE_ABRICATE as MERGE_ABRICATE_PLM; COMBINE_ABRICATE_RESULTS; PLOT_SUMMARIZE_ABRICATE; MLST_CONTIGS_PS; MERGE_MLST; BAKTA; PANAROO; RAXML_NG; IQTREE2; MULTIQC; KRAKEN2; GTDBTK_CLASSIFY } from './modules/nanopore.nf'

// -------- Params --------
params.kraken2_db = params.kraken2_db ?: null
params.gtdbtk_db  = params.gtdbtk_db  ?: null
params.reads  = params.reads ?: null
params.outdir = params.outdir ?: 'results'

// -------- Channels --------

def P = (params.reads ?: '').toString()
log.info "CWD: ${java.nio.file.Paths.get('').toAbsolutePath()}"
log.info "reads pattern: ${P}"

// JVM-side checks (no globbing, just: can Java see this path?)
def jf = new java.io.File(P)
log.info "JavaFile.exists=${jf.exists()} isFile=${jf.isFile()} path=${jf.getPath()}"

import java.nio.file.*
try {
  def np = Paths.get(P)
  log.info "NIO Files.exists=${Files.exists(np)} isRegularFile=${Files.isRegularFile(np)}"
} catch(Exception e) {
  log.info "NIO check threw: ${e.class.name}: ${e.message}"
}
if( params.reads ) {
  Channel.fromPath(params.reads, checkIfExists: true) // (sample_id, fastq.gz)
         .ifEmpty { error "No reads found for: ${params.reads}" }
         .map { f ->
         def sid = f.name
              .replaceFirst(/\.fastq\.gz$/, '')
              .replaceFirst(/\.fq\.gz$/, '')
              .replaceFirst(/\.fastq$/, '')
              .replaceFirst(/\.fq$/, '')
              tuple(sid, f)
              } 
         .set { reads } 
}
// Optional: drop empty files early
// reads = reads.filter { sid, fq -> fq.size() > 0 }

if( params.assemblies ) {
  Channel
    .fromPath(params.assemblies, checkIfExists: true)   // (sample_id, fasta)
    .ifEmpty { exit 1, "No assemblies found for: ${params.assemblies}" }
    .map { f ->
      def sid = f.getBaseName()
        .replaceFirst(/\.fa(sta)?$/, '')
        .replaceFirst(/\.fna$/, '')
      tuple(sid, f)
    }
    .set { assemblies }
}


// ---------------- Workflows ----------------

// -------- Workflow: assembly_amr_pangenome --------
// Takes nanopore reads (tuple: sample_id, reads), does QC, filtering, assembly,
// assembly evaluation, AMR and plasmid gene detection with ABRicate, annotation with Bakta,
// pangenome analysis with Panaroo, and generates a MultiQC report.

workflow assembly_amr_pangenome {

  main:
    // QC -> Filter -> QC
    np_raw   = NANOPLOT_RAW(reads)
    filt     = FILTLONG(reads)            // -> (sample_id, .filt.fastq.gz)
    np_filt  = NANOPLOT_FILT(filt.reads)

    //Taxonomy from filtered reads (contamination check)
    //kraken = KRAKEN2(filt.reads)

    // Assemble
    flye     = FLYE(filt.reads)           // -> (sample_id, .contigs.fasta)

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
    abri_summ_plots = PLOT_SUMMARIZE_ABRICATE(combined)

    // MLST (per sample, merge then summarise)
    mlst = MLST_CONTIGS_PS(flye.assembly, params.mlst_scheme)
    mlst_merged  = mlst.mlst_tsv.collect()
    mlst_summary = mlst.mlst_summary(mlst_merged) 

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
                //.mix( kraken.report.map { _, f -> f } )
                .mix( tag_amr.tagged )
                .mix( tag_plm.tagged )
                .collect()
    mqc = MULTIQC(mqc_inputs)
    
  emit:
    nanoplot_raw_dirs    = np_raw.report_dir
    nanoplot_filt_dirs   = np_filt.report_dir
    filtered_reads       = filt.reads
    //kraken_report        = kraken.report
    //kraken_calls         = kraken.classified
    assemblies           = flye.assembly
    gfa_files            = flye.graph
    asm_info_files       = flye.info
    quast_dir            = quast.report_dir
    quast_tsv            = quast.report_tsv
    quast_txt            = quast.report_txt
    bandage_png         = bandage.png
    bandage_svg         = bandage.svg
    bandage_nodes_table = bandage.info
    abricate_amr_per_sample  = ab.amr_tsv
    abricate_amr_summary     = ab.amr_summary
    abricate_plm_per_sample  = ab.plm_tsv
    abricate_plm_summary     = ab.plm_summary
    abricate_all_per_sample  = abri_all
    abricate_merged_tsv      = abri_merged
    abricate_combined        = combined.combined
    abricate_summary_plots   = abri_summ_plots.plots
    per_sample_summary       = abri_summ_plots.per_sample_summary
    mlst_tsv                 = mlst.mlst_tsv
    mlst_merged              = mlst.mlst_merged
    mlst_summary             = mlst.mlst_summary
    //bakta_gff            = bakta.gff
    //core_alignment       = pana.core_alignment
    //raxml_tree           = raxml.besttree
    //iqtree_tree          = iq.treefile
    multiqc_report       = mqc.report
}


// -------- Workflow: taxonomy_from_reads --------
// Takes nanopore reads (tuple: sample_id, reads), runs Kraken2 on reads,
// assembles with Flye, then runs GTDB-Tk on the assembly.
workflow taxonomy_from_reads {

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
// workflow { taxonomy_from_reads() }

// -------- Workflow: amr_pangenome_from_assemblies --------
// Takes assemblies (tuple: sample_id, assembly), runs ABRicate for AMR genes, plasmid genes,
// annotates with Bakta, does pangenome analysis with Panaroo, and builds trees
// from the core alignment with RAxML-NG and IQ-TREE2.

workflow amr_pangenome_from_assemblies {

  main:
    asm = assemblies

    // Assembly evaluation
    quast = QUAST(asm)

    // Summarise Quast results with MultiQC
    mqc = MULTIQC(quast.report_dir.map { _, d -> d }.collect())
    

  // ABRicate (per sample), tag, merge
    ab          = ABRICATE(asm)
    tag_amr     = ABRICATE_TAG_AMR(ab.amr_tsv, 'amr')
    tag_plm     = ABRICATE_TAG_PLM(ab.plm_tsv, 'plm')
    merged_amr  = MERGE_ABRICATE_AMR(tag_amr.tagged.collect(), 'amr')
    merged_plm  = MERGE_ABRICATE_PLM(tag_plm.tagged.collect(), 'plm')
    combined    = COMBINE_ABRICATE_RESULTS(merged_amr.merged, merged_plm.merged)
    abri_all    = ab.amr_tsv.mix(ab.plm_tsv).collect()
    abri_merged = merged_amr.merged.mix(merged_plm.merged).collect()

  // Bakta (per sample)
  //bakta    = BAKTA_ANNOTATE(asm)

  // ----- Collect GFFs for cohort-level steps -----
  //gff_list = bakta.gff
              //.map { sid, gff -> gff }
              //.collect()
              .map { gffs -> tuple(asm_ch.map{sid,_->sid}.toList().unique(), gffs) }

  // Panaroo + Piggy on the whole cohort
  //panaroo = PANAROO_RUN(gff_list)
 //piggy   = PIGGY_RUN(gff_list)

  // ----- Trees from Panaroo core alignment -----
  // Run both trees only if a core alignment was produced
  //raxml  = RAXML_NG_TREE(panaroo.core_aln)
  //iqtree = IQTREE2_TREE(panaroo.core_aln)

  // ----- Emits -----
  emit:
    quast_dir            = quast.report_dir
    quast_tsv            = quast.report_tsv
    quast_txt            = quast.report_txt
    multiqc_report       = mqc.report
    abricate_amr_per_sample  = ab.amr_tsv
    abricate_amr_summary     = ab.amr_summary
    abricate_plm_per_sample  = ab.plm_tsv
    abricate_plm_summary     = ab.plm_summary
    abricate_all_per_sample  = abri_all
    abricate_merged_tsv  = abri_merged
    abricate_combined    = combined.combined
    //gff_files    = bakta.gff
    //panaroo_dir  = panaroo.outdir
    //piggy_dir    = piggy.outdir
    //raxml_files  = raxml.files
    //iqtree_files = iqtree.files
}