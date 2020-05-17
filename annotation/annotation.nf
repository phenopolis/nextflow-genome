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

Vep_ch = For_vep_ch.join(Cadd_ch)

vep_human_ref_bundle = "${human_ref}.gz ${human_ref}.fai"

PLUGIN_NAMES = params.vep_plugins.collect { it -> it.name }
PLUGIN_ENVS = params.vep_plugins.collect { it -> it.env }
PLUGIN_CLOUDS = params.vep_plugins.collect { it -> it.cloud }
PLUGIN_LOCALS = params.vep_plugins.collect { it -> it.local }
PLUGIN_EXTRAS = params.vep_plugins.collect { it -> it.extra }
PLUGIN_FORMATS = params.vep_plugins.collect { it -> it.format }

ANNOTATION_ENVS = params.vep_annotations.collect { it -> it.env }
ANNOTATION_CLOUDS = params.vep_annotations.collect { it -> it.cloud }
ANNOTATION_LOCALS = params.vep_annotations.collect { it -> it.local }
ANNOTATION_ANNOTATIONS = params.vep_annotations.collect { it -> it.annotations }
ANNOTATION_FORMATS = params.vep_annotations.collect { it -> it.format }

process 'vep' {
  tag "${targetId}"
  //publishDir "."
  cpus params.vep_threads
  label 'small_batch'
  memory '15 G'
  //errorStrategy 'retry'
  //maxRetries 2
  container params.vep_docker
  echo true
  maxForks 1

  input:
    path(cadd_file) from Cadd_ch
    path("*") from CaddIndex_ch
    tuple val(targetId), path(input_vcf), path(input_bed) from For_vep_ch
  output:
    path "*.vep.json" into Vep_ch

  when:
    targetId == 'exome_target_00'

  script:
    def plugin_names = '"' + PLUGIN_NAMES.join('" "') + '"'
    def plugin_envs = '"' + PLUGIN_ENVS.join('" "') + '"'
    def plugin_clouds = '"' + PLUGIN_CLOUDS.join('" "') + '"'
    def plugin_locals = '"' + PLUGIN_LOCALS.join('" "') + '"'
    def plugin_extras = '"' + PLUGIN_EXTRAS.join('" "') + '"'
    def plugin_formats = '"' + PLUGIN_FORMATS.join('" "') + '"'
    def annotation_envs = '"' + ANNOTATION_ENVS.join('" "') + '"'
    def annotation_clouds = '"' + ANNOTATION_CLOUDS.join('" "') + '"'
    def annotation_locals = '"' + ANNOTATION_LOCALS.join('" "') + '"'
    def annotation_annotations = '"' + ANNOTATION_ANNOTATIONS.join('" "') + '"'
    def annotation_formats = '"' + ANNOTATION_FORMATS.join('" "') + '"'

    """
    source s3.bash
    source cadd_download_helper.sh
    aws_profile="${params.aws_ref_profile}"

    plugin_names=(${plugin_names})
    plugin_envs=(${plugin_envs})
    plugin_clouds=(${plugin_clouds})
    plugin_locals=(${plugin_locals})
    plugin_extras=(${plugin_extras})
    plugin_formats=(${plugin_formats})
    annotation_envs=(${annotation_envs})
    annotation_clouds=(${annotation_clouds})
    annotation_locals=(${annotation_locals})
    annotation_annotations=(${annotation_annotations})
    annotation_formats=(${annotation_formats})

    chrom=\$(head -1 ${input_bed} | cut -f1)
    start=\$(head -1 ${input_bed} | cut -f2)
    start=\$(( \${start}-${params.input_padding} ))
    start=\$(( \${start} > 0 ? \${start} : 0 ))
    end=\$(tail -1 ${input_bed} | cut -f3)
    end=\$(( \${end}+${params.input_padding} ))

    downloads=()
    # human_refs
    for downfile in ${vep_human_ref_bundle}
    do
      downloads+=("nxf_s3_retry nxf_s3_download ${human_ref_path}/\$downfile ./\$downfile")
    done
    # download annotation databases
    for ind in \${!annotation_clouds[@]}
    do
      cloud=\$(printf "\${annotation_clouds[\$ind]}" "\$chrom")
      local=\${annotation_locals[\$ind]}
      mkdir -p \$(dirname \${local})
      ls \$(dirname \${local})
      if [[ "\$cloud" == *gz ]]
      then
        downloads+=("\${annotation_envs[\$ind]} tabix -h \$cloud \${chrom}:\${start}-\${end} | bgzip -c > \$local && tabix \${annotation_formats[\$ind]} \$local")
      else
        downloads+=("nxf_s3_retry nxf_s3_download \$cloud \$local")
      fi
    done
    # download plugin script
    mkdir -p ~/.vep/Plugins
    aws_profile="${params.vep_s3_profile}"
    downloads+=("nxf_s3_retry nxf_s3_download ${params.vep_plugins_path} ~/.vep/Plugins")

    # download cache
    # download relevant files in vep
    vep_files=\$(/home/ec2-user/miniconda/bin/aws s3 ls ${params.vep_s3_profile} ${params.vep_database_path}/${params.vep_cache}/\${chrom}/ | awk '{print \$4'})
    for vep_file in \${vep_files[@]}
    do
      is_overlap=\$(overlap \$start \$end \${vep_file})
      if [ \${is_overlap} -eq 1 ]; then
          source="${params.vep_database_path}/${params.vep_cache}/\${chrom}/\$vep_file"
          target="~/.vep/${params.vep_cache}/\${chrom}/\$vep_file"
          target_dir=\$(dirname "\${target}")
          mkdir -p \${target_dir}
          downloads+=("nxf_s3_retry nxf_s3_download \${source} \${target}")
      fi
    done
    
    # download rest of vep files not in the vep_chrom list
    vep_chroms=(${params.vep_chromosomes})
    pattern=\$(printf 'PRE %s/\\|' "\${vep_chroms[@]}") 
    pattern=\${pattern%\\\\|}
    inclusion=(\$(/home/ec2-user/miniconda/bin/aws s3 ls ${params.vep_s3_profile} ${params.vep_database_path}/${params.vep_cache}/ | grep -v "\$pattern" | awk '{print \$NF}'))
    inclusion=(\$(printf "${params.vep_cache}/%s " "\${inclusion[@]}"))
    mkdir -p "~/.vep/${params.vep_cache}"
    for include in \${inclusion[@]}
    do
      downloads+=("nxf_s3_retry nxf_s3_download ${params.vep_database_path}/\${include} ~/.vep/\${include}")
    done

    # download plugin databases
    for ind in \${!plugin_clouds[@]}
    do
      cloud=\${plugin_clouds[\$ind]}
      local=\${plugin_locals[\$ind]}
      if [[ "\$cloud" == *gz ]]
      then
        mkdir -p \$(dirname \${local})
        comm="\${plugin_envs[\$ind]} tabix -h \$cloud \${chrom}:\${start}-\${end} | bgzip -c > \$local && tabix \${plugin_formats[\$ind]} \$local"
        downloads+=("\$comm")
      elif [ ! -z "\$cloud" ]
      then
        mkdir -p \$(dirname \${local})
        downloads+=("nxf_s3_retry nxf_s3_download \$cloud \$local")
      fi
    done

    # shuffle downloads since tabix will struggle
    declare -a outarray
    perm "\${downloads[@]}"

    # a plugin uses Math::CDF that is not in place
    downloads+=("cpanm Math::CDF")

    nxf_parallel "\${outarray[@]}"

    # unzip human_ref
    gunzip ${human_ref}.gz

    # vep
    plugins=
    for ind in \${!plugin_clouds[@]}
    do
      plugin_local=
      if [ ! -z "\${plugin_locals[\$ind]}" ]
      then
        plugin_local=",\${plugin_locals[\$ind]}"
      fi
      plugins="\${plugins} --plugin \${plugin_names[\$ind]}\${plugin_local}\${plugin_extras[\$ind]}"
    done
    # add cadd
    plugins="\${plugins} --plugin CADD,${cadd_file}"

    annotations=
    for ind in \${!annotation_clouds[@]}
    do
      localfile=\${annotation_locals[\$ind]}
      any_entry=\$(zcat \${localfile} | grep -m 1 -v # -c)
      [[ \${any_entry} == 1 ]] && annotations="\${annotations} --custom \${localfile},\${annotation_annotations[\$ind]}"
    done
    annotations="\${annotations} --custom ${input_vcf},input,vcf,exact,0,${params.SELF_INFO_FIELDS}"

    vep -i ${input_vcf} ${params.vep_flags} --fork ${params.vep_threads} --dir ~/.vep --fasta ${human_ref} -a ${params.build} \$annotations \$plugins -o a.vep.json

    """

}
