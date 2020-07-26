HG38 = ['hg38', 'GRCh38', 'grch38']
HG19 = ['hg19', 'GRCh37', 'b37', 'grch37']
human_ref_base = null
human_ref_path = null
if (HG38.contains(params.build)) {
  human_ref_base = params.human_ref_base_hg38
  human_ref_path = params.human_ref_path_hg38
} else if (HG19.contains(params.build)) {
  human_ref_path = params.human_ref_path_hg19
  human_ref_base = params.human_ref_base_hg19
} else {
  exit 1, "cannot understand build. Please provide either --build hg19, or --build hg38 (default is hg19)"
}
human_ref = "${human_ref_base}.fasta"


Channel.fromList(file(params.input_table).readLines().collect{[it, "${params.input_path}/$it/${params.input_filename}", "${params.input_path}/$it/${params.input_filename}.tbi" ] }).into {Cadd_input_ch; Vep_input_ch}

process 'cadd' {
  tag "${sampleId}"
  cpus 1
  label 'cadd'
  memory '7 G'

  container params.cadd_docker

  input:
    tuple val(sampleId), path(input_vcf), path(input_vcf_tbi) from Cadd_input_ch
  output:
    tuple sampleId, path("cadd.tsv.gz"), path("cadd.tsv.gz.tbi") into Cadd_ch

  """

  # run cadd
  cadd ${params.cadd_flags} -o cadd.tsv.gz ${input_vcf} && tabix -p vcf cadd.tsv.gz
  """
}

Vep_in_ch = Vep_input_ch.join(Cadd_ch)

vep_human_ref_bundle = "${human_ref}.gz ${human_ref}.fai"

PLUGIN_NAMES = params.vep_plugins.collect { it -> it.name }
PLUGIN_LOCALS = params.vep_plugins.collect { it -> it.local }
PLUGIN_EXTRAS = params.vep_plugins.collect { it -> it.extra }

ANNOTATION_LOCALS = params.vep_annotations.collect { it -> it.local }
ANNOTATION_ANNOTATIONS = params.vep_annotations.collect { it -> it.annotations }

process 'vep' {
  tag "${sampleId}"
  cpus params.vep_threads
  label 'vep'
  memory '25 G'
  //errorStrategy 'retry'
  //maxRetries 2
  container params.vep_docker

  input:
    tuple val(sampleId), path(input_vcf), path(input_vcf_index), path(input_cadd), path(cadd_index) from Vep_in_ch
  output:
    path "vep.json" into Vep_ch

  script:
    def plugin_names = '"' + PLUGIN_NAMES.join('" "') + '"'
    def plugin_locals = '"' + PLUGIN_LOCALS.join('" "') + '"'
    def plugin_extras = '"' + PLUGIN_EXTRAS.join('" "') + '"'
    def annotation_locals = '"' + ANNOTATION_LOCALS.join('" "') + '"'
    def annotation_annotations = '"' + ANNOTATION_ANNOTATIONS.join('" "') + '"'


    """
    source s3.bash

    plugin_names=(${plugin_names})
    plugin_locals=(${plugin_locals})
    plugin_extras=(${plugin_extras})
    annotation_locals=(${annotation_locals})
    annotation_annotations=(${annotation_annotations})

    aws_profile="${params.aws_ref_profile}"

    downloads=()
    # human_refs
    for downfile in ${vep_human_ref_bundle}
    do
      downloads+=("nxf_s3_retry nxf_s3_download ${human_ref_path}/\$downfile ./\$downfile")
    done

    nxf_parallel "\${downloads[@]}"

    # a plugin uses Math::CDF that is not in place
    cpanm Math::CDF

    # unzip human_ref
    gunzip ${human_ref}.gz

    # vep
    plugins=
    for ind in \${!plugin_names[@]}
    do
      plugin_local=
      if [ ! -z "\${plugin_locals[\$ind]}" ]
      then
        plugin_local=",\${plugin_locals[\$ind]}"
      fi
      plugins="\${plugins} --plugin \${plugin_names[\$ind]}\${plugin_local}\${plugin_extras[\$ind]}"
    done
    # add cadd
    plugins="\${plugins} --plugin CADD,${input_cadd}"

    annotations=
    for ind in \${!annotation_locals[@]}
    do
      localfile=\${annotation_locals[\$ind]}
      annotations="\${annotations} --custom \${localfile},\${annotation_annotations[\$ind]}"
    done
    annotations="\${annotations} --custom ${input_vcf},input,vcf,exact,0,${params.SELF_INFO_FIELDS}"

    vep -i ${input_vcf} ${params.vep_flags} --fork ${params.vep_threads} --dir /data/.vep --dir_cache /data/.vep --fasta ${human_ref} -a ${params.build} \$annotations \$plugins -o vep.json

    uploads=("nxf_s3_retry nxf_s3_upload vep.json ${params.s3_deposit}/${sampleId}")
    aws_profile="${params.s3_deposit_profile}"
    nxf_parallel "\${uploads[@]}"


    """

}

