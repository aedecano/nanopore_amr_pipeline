nextflow.enable.dsl = 2

// ------------- NanoPlot (raw reads) -------------
process NANOPLOT_RAW {
  tag "$sample_id"
  label 'light'
  //conda 'bioconda::nanoplot=1.42.0'
  //conda 'conda_setup/envs/qc.yaml'
  conda "$HOME/.conda_envs_nf_cache/env-nf-qc"
  publishDir "${params.outdir}/nanoplot_raw", mode: 'copy'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}_nanoplot_raw"), emit: report_dir

  script:
  """
  set -eo pipefail
  mkdir -p "${sample_id}_nanoplot_raw"
  NanoPlot --fastq "${reads}" --N50 --loglength --threads ${task.cpus} -o "${sample_id}_nanoplot_raw"
  """
}

// ------------- Filtlong -------------
process FILTLONG {
  tag "$sample_id"
  label 'light'
  //conda 'bioconda::filtlong=0.2.1'
  conda "$HOME/.conda_envs_nf_cache/env-filtlong"
  publishDir "${params.outdir}/filtlong", mode: 'copy'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}.filt.fastq.gz"), emit: reads

  script:
  """
  set -eo pipefail
  test -s "${reads}" || { echo "Input reads missing/empty: ${reads}" >&2; exit 1; }

  filtlong --min_length 1000 --target_bases 200000000 "${reads}" \
    | gzip -c > "${sample_id}.filt.fastq.gz"
  """
}

// ------------- NanoPlot (post-Filtlong) -------------
process NANOPLOT_FILT {
  tag "$sample_id"
  label 'light'
  //conda 'bioconda::nanoplot=1.42.0'
  //conda 'conda_setup/envs/qc.yaml'
  conda "$HOME/.conda_envs_nf_cache/env-nf-qc"
  publishDir "${params.outdir}/nanoplot_filt", mode: 'copy'

  input:
    tuple val(sample_id), path(filt)

  output:
    tuple val(sample_id), path("${sample_id}_nanoplot_filt"), emit: report_dir

  script:
  """
  set -eo pipefail
  mkdir -p "${sample_id}_nanoplot_filt"
  NanoPlot --fastq "${filt}" --N50 --loglength --threads ${task.cpus} -o "${sample_id}_nanoplot_filt"
  """
}

// ------------- Flye assembly -------------
process FLYE {
  tag "$sample_id"
  label 'medium'
  conda "$HOME/.conda_envs_nf_cache/env-flye"
  publishDir "${params.outdir}/flye/${sample_id}", mode: 'copy'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}.contigs.fasta"),        emit: assembly
    tuple val(sample_id), path("${sample_id}_graph.gfa"),  optional: true, emit: graph
    tuple val(sample_id), path("${sample_id}.assembly_info.txt"), optional: true, emit: info

  script:
  """
  set -euo pipefail

  test -s "${reads}" || { echo "Filtered reads missing/empty: ${reads}" >&2; exit 1; }

  flye --nano-raw "${reads}" --out-dir . --threads ${task.cpus}

  # Standardize filenames for downstream rules
  mv assembly.fasta "${sample_id}.contigs.fasta"

  if [[ -s assembly_graph.gfa ]]; then
    cp assembly_graph.gfa "${sample_id}_graph.gfa"
  fi

  if [[ -s assembly_info.txt ]]; then
    cp assembly_info.txt "${sample_id}.assembly_info.txt"
  fi
  """
}

// ------------- QUAST (assembly evaluation) -------------
process QUAST {
  tag "$sample_id"
  label 'medium'
  
  conda "$HOME/.conda_envs_nf_cache/env-quast-x64"
  publishDir "${params.outdir}/quast/${sample_id}", mode: 'copy'

  input:
    tuple val(sample_id), path(assembly)

  output:
    tuple val(sample_id), path("${sample_id}_quast"),                            emit: report_dir
    tuple val(sample_id), path("${sample_id}_quast/report.tsv"), optional: true, emit: report_tsv
    tuple val(sample_id), path("${sample_id}_quast/report.txt"), optional: true, emit: report_txt

  script:
  """
  set -euo pipefail
  mkdir -p ${sample_id}_quast

  #if command -v quast.py >/dev/null 2>&1; then
    #quast.py "${assembly}" -o "${sample_id}_quast" -t ${task.cpus} --no-check --silent
  #else
    quast    "${assembly}" -o "${sample_id}_quast" -t ${task.cpus} --no-check --silent
  #fi
  """
}

// ------------- Bandage (export images & node table) -------------
process BANDAGE_IMAGE {
  tag "$sample_id"
  label 'graph'
  publishDir "${params.outdir}/bandage/${sample_id}", mode: 'copy'

  input:
    tuple val(sample_id), path(gfa)

  output:
    tuple val(sample_id), path("${sample_id}_graph.png"),                 emit: png
    tuple val(sample_id), path("${sample_id}_graph.svg"),  optional: true, emit: svg
    tuple val(sample_id), path("${sample_id}_graph_nodes.txt"), optional: true, emit: info

  script:
  """
  set -euo pipefail
  test -s "${gfa}" || { echo "Missing or empty GFA: ${gfa}" >&2; exit 1; }

  # PNG export
  Bandage image "${gfa}" "${sample_id}_graph.png" --resolution ${params.bandage_resolution:-300} ${params.bandage_annotations ? '--annotations' : ''}

  # SVG export (skip if build lacks --format)
  if Bandage image --help 2>/dev/null | grep -q -- '--format'; then
    Bandage image "${gfa}" "${sample_id}_graph.svg" --format svg ${params.bandage_annotations ? '--annotations' : ''}
  fi

  # Node/coverage table export (optional)
  Bandage info "${gfa}" "${sample_id}_graph_nodes.txt" || true
  """
}

// ------------- ABRicate (per sample) -------------
process ABRICATE {
  tag "$sample_id"
  label 'light'
  //conda 'bioconda::abricate=1.0.1'
  publishDir "${params.outdir}/abricate", mode: 'copy'

  input:
    tuple val(sample_id), path(contigs)

  output:
    tuple val(sample_id), path("${sample_id}.abricate.tsv"),         emit: tsv
    tuple val(sample_id), path("${sample_id}.abricate.summary.tsv"), emit: summary

  when:
    params.abricate_db != null

  script:
  """
  set -eo pipefail
  abricate --db ${params.abricate_db} "${contigs}" > "${sample_id}.abricate.tsv"
  abricate --summary "${sample_id}.abricate.tsv" > "${sample_id}.abricate.summary.tsv"
  """
}

// ------------- ABRicate tag sample_id for merge -------------
process ABRICATE_TAG {
  tag "$sample_id"
  label 'light'
  publishDir "${params.outdir}/abricate", mode: 'copy'

  input:
    tuple val(sample_id), path(tsv)

  output:
    path("${sample_id}.abricate.tagged.tsv"), emit: tagged

  script:
  """
  set -eo pipefail
  # Write header with a new leading 'sample_id' column
  { printf "sample_id\\t"; head -n1 "${tsv}"; } > "${sample_id}.abricate.tagged.tsv"

  # Append body rows, prefixing sample_id to each line
  tail -n +2 "${tsv}" | sed "s/^/${sample_id}\\t/" >> "${sample_id}.abricate.tagged.tsv"
  """
}

// ------------- Merge all ABRicate tagged TSVs -------------
process MERGE_ABRICATE {
  tag "merge-abricate"
  label 'light'
  publishDir "${params.outdir}/abricate", mode: 'copy'

  input:
    path(tagged_tsvs)

  output:
    path("abricate_merged.tsv"), emit: merged

  script:
  """
  set -eo pipefail
  first=\$(ls -1 ${tagged_tsvs} | head -n1)
  { head -n1 "\$first"; for f in ${tagged_tsvs}; do awk 'NR>1' "\$f"; done; } > abricate_merged.tsv
  """
}

// ------------- Bakta (functional annotation) -------------
process BAKTA {
  tag "$sample_id"
  label 'medium'
  conda 'bioconda::bakta=1.9.4'
  publishDir "${params.outdir}/bakta", mode: 'copy'

  input:
    tuple val(sample_id), path(contigs)

  output:
    tuple val(sample_id), path("${sample_id}.gff3"), emit: gff
    tuple val(sample_id), path("${sample_id}.gbff"), emit: gbff
    tuple val(sample_id), path("${sample_id}.faa"),  emit: proteins
    tuple val(sample_id), path("${sample_id}.ffn"),  emit: genes

  when:
    params.bakta_db != null

  script:
  """
  set -eo pipefail
  bakta "${contigs}" --db ${params.bakta_db} --prefix "${sample_id}" --output "${PWD}" --threads ${task.cpus}
  test -f "${sample_id}.gff" && mv "${sample_id}.gff" "${sample_id}.gff3" || true
  """
}

// ------------- Panaroo (pangenome) -------------
process PANAROO {
  tag "panaroo"
  label 'heavy'
  publishDir "${params.outdir}/panaroo", mode: 'copy'

  input:
    path(gff_list)

  output:
    path("panaroo_output/core_gene_alignment.aln.fasta"), emit: core_alignment
    path("panaroo_output/roary_output.csv"),               emit: roary_like
    path("panaroo_output/graph.gml"),                      emit: graph

  script:
  """
  set -eo pipefail
  mkdir -p panaroo_output
  panaroo -i ${gff_list.join(' ')} -o panaroo_output -t ${task.cpus} --clean-mode strict --remove-invalid-genes
  """
}

// ------------- RAxML-NG (tree) -------------
process RAXML_NG {
  tag "raxml"
  label 'heavy'
  conda 'bioconda::raxml-ng=1.2.2'
  publishDir "${params.outdir}/phylo", mode: 'copy'

  input:
    path(aln)

  output:
    path("raxml.bestTree"),  emit: besttree
    path("raxml.bestModel"), emit: bestmodel
    path("raxml.log"),       emit: log

  script:
  """
  set -eo pipefail
  raxml-ng --msa "${aln}" --all --model GTR+G --seed 13 --bs-trees 100 --threads ${task.cpus} --prefix raxml
  """
}

// ------------- IQ-TREE2 (alternate tree) -------------
process IQTREE2 {
  tag "iqtree2"
  label 'heavy'
  conda 'bioconda::iqtree=2.3.5'
  publishDir "${params.outdir}/phylo_iqtree", mode: 'copy'

  input:
    path(aln)

  output:
    path("core.aln.fasta.treefile"), emit: treefile
    path("core.aln.fasta.iqtree"),   emit: iqreport
    path("core.aln.fasta.log"),      emit: log

  script:
  """
  set -eo pipefail
  cp "${aln}" core.aln.fasta
  iqtree2 -s core.aln.fasta -m MFP -bb 1000 -alrt 1000 -nt ${task.cpus} -redo
  """
}

// ------------- MultiQC (aggregate) -------------
process MULTIQC {
  tag "multiqc"
  label 'light'
  //conda 'bioconda::multiqc=1.25'
  publishDir "${params.outdir}/multiqc", mode: 'copy'

  input:
    path(report_inputs)

  output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"),        emit: data_dir

  script:
  """
  set -eo pipefail
  mkdir -p mqc_out
  multiqc -o mqc_out ${report_inputs.join(' ')}
  mv mqc_out/multiqc_report.html .
  mv mqc_out/multiqc_data .
  """
}

// ------------- Kraken2 (taxonomy on reads) -------------
process KRAKEN2 {
  tag "$sample_id"
  label 'light'
  
  //conda 'bioconda::kraken2=2.1.3'
  conda '$HOME/miniforge3/envs/kraken2'
  publishDir "${params.outdir}/kraken2", mode: 'copy'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}.kraken2.report"), emit: report
    tuple val(sample_id), path("${sample_id}.kraken2.tsv"),    emit: classified
    // If you also produce a Krona HTML, add it as a separate line:
    // tuple val(sample_id), path("${sample_id}.krona.html"),    emit: krona

  when:
    params.kraken2_db != null

  script:
  """
  set -eo pipefail
  kraken2 --db ${params.kraken2_db} \
          --threads ${task.cpus} \
          --report "${sample_id}.kraken2.report" \
          --output "${sample_id}.kraken2.tsv" \
          "${reads}"

  # Optional: Krona (if installed and desired)
  # cut -f2,3 "${sample_id}.kraken2.tsv" | ktImportTaxonomy -o "${sample_id}.krona.html" -
  """
}

// ------------- GTDB-Tk (taxonomy on assemblies) -------------
process GTDBTK_CLASSIFY {
  tag "$sample_id"
  label 'heavy'
  conda 'bioconda::gtdbtk=2.4.0'
  publishDir "${params.outdir}/gtdbtk", mode: 'copy'

  input:
    tuple val(sample_id), path(assembly)

  output:
    tuple val(sample_id), path("${sample_id}_gtdbtk_summary.tsv"),      emit: summary
    tuple val(sample_id), path("${sample_id}_gtdbtk_classification.tsv"), emit: classification
    tuple val(sample_id), path("${sample_id}_gtdbtk_logs"),               emit: logs_dir

  when:
    params.gtdbtk_db && params.gtdbtk_db != ""

  script:
  """
  set -eo pipefail
  mkdir -p ${sample_id}_gtdbtk_work
  gtdbtk classify_wf \\
    --genome ${assembly} \\
    --out_dir ${sample_id}_gtdbtk_work \\
    --data_dir ${params.gtdbtk_db} \\
    --cpus ${task.cpus ?: 8}
  # Collect key outputs with consistent names
  cp ${sample_id}_gtdbtk_work/classify/summary.tsv ${sample_id}_gtdbtk_summary.tsv || true
  cp ${sample_id}_gtdbtk_work/classify/gtdbtk.bac120.classification.tsv ${sample_id}_gtdbtk_classification.tsv || true
  cp -r ${sample_id}_gtdbtk_work/logs ${sample_id}_gtdbtk_logs || true
  """
}
