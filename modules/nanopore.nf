nextflow.enable.dsl = 2

// ------------- NanoPlot (raw reads) -------------
process NANOPLOT_RAW {
  tag "$sample_id"
  label 'light'
  publishDir "${params.outdir}/nanoplot_raw", mode: 'copy'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}_nanoplot_raw"), emit: report_dir

  script:
  """
  set -eo pipefail
  mkdir -p "${sample_id}_nanoplot_raw"
  NanoPlot --fastq "${reads}" --n50 --loglength --threads ${task.cpus} -o "${sample_id}_nanoplot_raw"
  """
}

// ------------- Filtlong -------------
process FILTLONG {
  tag "$sample_id"
  label 'light'
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
  publishDir "${params.outdir}/nanoplot_filt", mode: 'copy'

  input:
    tuple val(sample_id), path(filt)

  output:
    tuple val(sample_id), path("${sample_id}_nanoplot_filt"), emit: report_dir

  script:
  """
  set -eo pipefail
  mkdir -p "${sample_id}_nanoplot_filt"
  NanoPlot --fastq "${filt}" --n50 --loglength --threads ${task.cpus} -o "${sample_id}_nanoplot_filt"
  """
}

// ------------- Flye assembly -------------
process FLYE {
  tag "$sample_id"
  label 'heavy'
  publishDir "${params.outdir}/flye", mode: 'copy'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}.contigs.fasta"), emit: assembly

  script:
  """
  set -eo pipefail
  test -s "${reads}" || { echo "Filtered reads missing/empty: ${reads}" >&2; exit 1; }

  flye --nano-raw "${reads}" --out-dir . --threads ${task.cpus}
  mv assembly.fasta "${sample_id}.contigs.fasta"
  """
}

// ------------- ABRicate (per sample) -------------
process ABRICATE {
  tag "$sample_id"
  label 'light'
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