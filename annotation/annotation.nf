/*
nextflow run annotation.nf --input_table input_table.tsv -bucket-dir s3://phenopolis-nextflow/nextflow-work/annotation -with-report annotation.html [-resume]

CADD should take around 10 minutes to finish.
VEP takes around 50 minutes (16 threads). Increase vep_threads in nextflow.config to boost speed
*/
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


//Channel.fromPath(params.input_table)
//  .splitCsv(header:['sampleId', 'read1', 'read2'], sep:',')
//  .map{ row -> tuple(row.sampleId, "${params.input_path}/${row.sampleId}/${params.input_filename}",  "${params.input_path}/${row.sampleId}/${params.input_filename}.tbi") }
Channel.fromList(file(params.input_table).readLines().collect{ it -> it.tokenize(',')[0]})
  .set {Input_ch}

Bed_ch = Channel.value(file("${params.input_beds}/*.bed"))

process 'horizontal_split' {
  tag "${sampleId}"
  cpus 4
  label 'small_batch'
  memory '8 G'
  container params.gatk_docker
  input:
    val(sampleId) from Input_ch
    path("*") from Bed_ch
  output:
    path("${sampleId}.*.vcf.gz") into Horizontal_split_ch
    tuple sampleId, path(params.input_filename), path("${params.input_filename}.tbi") into Normalise_input_ch

  """
  source s3.bash

  aws_profile="${params.s3_deposit_profile}"
  downloads=()
  for filename in ${params.input_filename} ${params.input_filename}.tbi; do
      cloudFile=${params.s3_deposit}/${sampleId}/\$filename
      downloads+=("nxf_s3_retry nxf_s3_download \$cloudFile \$filename")
  done
  nxf_parallel "\${downloads[@]}"
  for bed_file in *.bed
  do
    targetName=\${bed_file%".bed"}
    chrom=\$(head -1 \${bed_file} | cut -f1)
    start=\$(head -1 \${bed_file} | cut -f2)
    end=\$(tail -1 \${bed_file} | cut -f3)
    tabix -h ${params.input_filename} \$chrom:\$start-\$end | bgzip -c > ${sampleId}.\${targetName}.vcf.gz
  done
  """
}

Cadd_input_ch = Horizontal_split_ch
  .flatten()


//cadd_files = file(params.cadd_file_list).readLines().join(' ')
process 'cadd' {
  tag "${sampleId},${targetId}"
  cpus 1
  label 'cadd'
  memory '7 G'

  container params.cadd_docker

  input:
    path(input_vcf) from Cadd_input_ch
  output:
    tuple sampleId, path("${sampleId}.${targetId}.cadd.tsv.gz") into Cadd_ch

  script:
    sampleId = input_vcf.name.tokenize('.')[0]
    targetId = input_vcf.name.tokenize('.')[1]
    """
    cadd ${params.cadd_flags} -o ${sampleId}.${targetId}.cadd.tsv.gz ${input_vcf} && tabix -p vcf ${sampleId}.${targetId}.cadd.tsv.gz
    """
}

Merge_cadd_input_ch = Cadd_ch.groupTuple()

process 'merge_cadd' {
  tag "${sampleId}"
  cpus 1
  label 'small_batch'
  memory '7 G'
  container params.vep_docker

  input:
    tuple val(sampleId), path("*") from Merge_cadd_input_ch

  output:
    tuple sampleId, path("cadd.tsv.gz"), path("cadd.tsv.gz.tbi") into Merged_cadd_output_ch

  """
  infiles=(\$(ls *.gz))
  zcat \${infiles[0]} | head -2 > cadd.tsv
  for filename in "\${infiles[@]}"; do
    zcat \$filename | grep -v '#' >> cadd.tsv
  done
  bgzip cadd.tsv
  tabix -p vcf cadd.tsv.gz

  """

}

/*
vep doesn't work well with non-normalised vcf
*/
vep_human_ref_bundle = "${human_ref}.gz ${human_ref}.fai"
process 'normalise' {
  tag "${sampleId}"
  cpus 1
  label 'SSD'
  memory '8 G'
  container params.bcftools_docker
  input:
    tuple val(sampleId), path(input_vcf), path(input_vcf_index) from Normalise_input_ch
  output:
    tuple sampleId, path("normalised.vcf.gz") into Vep_input_ch

  script:
  """
  source s3.bash

  aws_profile="${params.aws_ref_profile}"
  downloads=()
  # human_refs
  for downfile in ${vep_human_ref_bundle}
  do
    downloads+=("nxf_s3_retry nxf_s3_download ${human_ref_path}/\$downfile ./\$downfile")
  done
  nxf_parallel "\${downloads[@]}"
  gunzip ${human_ref}.gz

  bcftools norm -Ou -m -any ${input_vcf} | bcftools norm -Oz -f ${human_ref} -o normalised.vcf.gz
  """
}
Vep_in_ch = Vep_input_ch.join(Merged_cadd_output_ch)

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
    tuple val(sampleId), path(input_vcf), path(input_cadd), path(cadd_index) from Vep_in_ch
  output:
    tuple sampleId, path("vep.json") into Vep_ch

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
    # vcf input
    tabix -p vcf ${input_vcf}

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
    
    gzip -c vep.json > vep.json.gz
    uploads=("nxf_s3_retry nxf_s3_upload vep.json.gz ${params.s3_deposit}/${sampleId}")
    aws_profile="${params.s3_deposit_profile}"
    nxf_parallel "\${uploads[@]}"

    """
}
