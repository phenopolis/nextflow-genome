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
Input_vcf_ch = Channel.fromPath(params.input_vcf).map { file -> tuple(file.name.take(file.name.lastIndexOf('.vqsr.vcf.gz')), file ) }

process 'vertical_split' {
  // separate genotypes from locations
  // it is a basic python script, with zcat/bgzip/tabix. Being lazy for now to use gatk's docker as it has all of these
  cpus 2
  label 'small_batch'
  memory '32 G'
  //errorStrategy 'retry'
  //maxRetries 2
  container params.gatk_docker
  publishDir "${params.outdir}/normalised-vcf/${cohortId}", mode: "copy", pattern: "*.tab"
  input:
    tuple val(cohortId), path(input_vcf) from Input_vcf_ch
  output:
    tuple cohortId, path("*.tab") into Sample_ch
    tuple cohortId, path("*_for_VEP.vcf.gz*") into Vertical_split_ch

  """
  zcat ${input_vcf} | multiallele_to_single_vcf.py \
  --headers_for_vep CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT \
  --file_for_vep ${cohortId}_for_VEP.vcf \
  --calls_file ${cohortId}_geno.tab \
  --depth_file ${cohortId}_depth.tab \
  --rowname_file ${cohortId}_rowname.tab \
  --colname_file ${cohortId}_colname.tab
  bgzip -f -c ${cohortId}_for_VEP.vcf > ${cohortId}_for_VEP.vcf.gz
  tabix -f -p vcf ${cohortId}_for_VEP.vcf.gz
  rm ${cohortId}_for_VEP.vcf
  """
}

Bed_ch = Channel.value(file("${params.input_beds}/*.bed"))

process 'horizontal_split' {
  cpus 4
  label 'small_batch'
  memory '8 G'
  container params.gatk_docker
  input:
    tuple val(cohortId), path(input_vcf) from Vertical_split_ch
    path("*") from Bed_ch
  output:
    path("*.vcf.gz") into Horizontal_split_ch

  """
  for bed_file in *.bed
  do
    targetName=\${bed_file%".bed"}
    chrom=\$(head -1 \${bed_file} | cut -f1)
    start=\$(head -1 \${bed_file} | cut -f2)
    end=\$(tail -1 \${bed_file} | cut -f3)
    tabix -h ${input_vcf} \$chrom:\$start-\$end | bgzip -c > \${targetName}.vcf.gz
  done
  """
}

Horizontal_split_ch
  .flatten()
  .map { file -> tuple(file.name.take(file.name.lastIndexOf('.vcf.gz')), file ) }
  .map { it -> tuple(it[0], file(it[1]), file("${params.input_beds}/${it[0]}.bed"))}
  .into {For_cadd_ch; For_vep_ch;}

//cadd_files = file(params.cadd_file_list).readLines().join(' ')

process 'cadd' {
  tag "${targetId}"
  cpus 1
  label 'cadd'
  memory '7 G'

  container params.cadd_docker

  input:
    tuple val(targetId), path(input_vcf), path(input_bed) from For_cadd_ch
  output:
    tuple targetId, path("${targetId}.cadd.tsv.gz"), path("${targetId}.cadd.tsv.gz.tbi") into Cadd_ch

  """

  # run cadd
  cadd ${params.cadd_flags} -o ${targetId}.cadd.tsv.gz ${input_vcf} && tabix -p vcf ${targetId}.cadd.tsv.gz
  """
}

Vep_in_ch = For_vep_ch.join(Cadd_ch)

vep_human_ref_bundle = "${human_ref}.gz ${human_ref}.fai"

PLUGIN_NAMES = params.vep_plugins.collect { it -> it.name }
PLUGIN_LOCALS = params.vep_plugins.collect { it -> it.local }
PLUGIN_EXTRAS = params.vep_plugins.collect { it -> it.extra }

ANNOTATION_LOCALS = params.vep_annotations.collect { it -> it.local }
ANNOTATION_ANNOTATIONS = params.vep_annotations.collect { it -> it.annotations }

process 'vep' {
  tag "${targetId}"
  publishDir "${baseDir}/vep", mode: 'copy'
  //publishDir "."
  cpus params.vep_threads
  label 'vep'
  memory '25 G'
  //errorStrategy 'retry'
  //maxRetries 2
  container params.vep_docker

  input:
    tuple val(targetId), path(input_vcf), path(input_bed), path(input_cadd), path(cadd_index) from Vep_in_ch
  output:
    path "${targetId}.vep.json" into Vep_ch
  when:
    targetId == 'exome_target_00'

  script:
    def plugin_names = '"' + PLUGIN_NAMES.join('" "') + '"'
    def plugin_locals = '"' + PLUGIN_LOCALS.join('" "') + '"'
    def plugin_extras = '"' + PLUGIN_EXTRAS.join('" "') + '"'
    def annotation_locals = '"' + ANNOTATION_LOCALS.join('" "') + '"'
    def annotation_annotations = '"' + ANNOTATION_ANNOTATIONS.join('" "') + '"'


    """
    source s3.bash
    tabix -p vcf ${input_vcf}

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

    vep -i ${input_vcf} ${params.vep_flags} --fork ${params.vep_threads} --dir /data/.vep --dir_cache /data/.vep --fasta ${human_ref} -a ${params.build} \$annotations \$plugins -o ${targetId}.vep.json
    

    """

}

