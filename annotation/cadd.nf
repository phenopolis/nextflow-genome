// to run: nextflow cadd.nf -bucket-dir s3://phenopolis-nextflow/nextflow-work -with-report cadd.html
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

Bed_ch = Channel.value(file("${params.input_beds}/*.bed"))
Input_ch = Channel.fromPath('/SAN/vyplab/UCLex/mainset_April2020/vulliamy_all2_for_VEP.fixed.vcf.gz')
Input_index_ch = Channel.fromPath('/SAN/vyplab/UCLex/mainset_April2020/vulliamy_all2_for_VEP.fixed.vcf.gz.tbi')
cohortId = 'vulliamy_all2'

process 'horizontal_split' {
  cpus 4
  memory '8 G'
  container params.gatk_docker
  input:
    path(input_vcf) from Input_ch
    path(input_vcf_index) from Input_index_ch
    path("*") from Bed_ch
  output:
    path("exome_target_*.vcf.gz") into Horizontal_split_ch

  """
  for bed_file in *.bed
  do
    targetName=\${bed_file%".bed"}
    chrom=\$(head -1 \${bed_file} | cut -f1)
    start=\$(head -1 \${bed_file} | cut -f2)
    end=\$(tail -1 \${bed_file} | cut -f3)
    tabix -h ${input_vcf} \$chrom:\$start-\$end > temp.vcf
    bgzip -c temp.vcf > \${targetName}.vcf.gz
    rm temp.vcf
  done
  """
}
Horizontal_split_ch
  .flatten()
  .map { file -> tuple(file.name.take(file.name.lastIndexOf('.vcf.gz')), file ) }
  .map { it -> tuple(it[0], file(it[1]), file("${params.input_beds}/${it[0]}.bed"))}
  .into {For_cadd_ch; For_vep_ch;}

process 'cadd' {
  tag "${targetId}"
  cpus 1
  label 'cadd'
  memory '7 G'
  publishDir "${baseDir}/cadds", mode: 'copy'

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
  publishDir "${baseDir}/vep", mode: 'copy'
  //publishDir "."
  cpus params.vep_threads
  label 'small_batch'
  memory '15 G'
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
    source utils.sh
    aws_profile="${params.aws_ref_profile}"
    tabix -p vcf ${input_vcf}

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

    echo "where is gsutil"
    ls /root/google-cloud-sdk/bin/gsutil
    downloads=()
    mkdir -p /data/.vep

    # human_refs
    for downfile in ${vep_human_ref_bundle}
    do
      downloads+=("nxf_s3_retry nxf_s3_download ${human_ref_path}/\$downfile ./\$downfile")
    done

    # download plugin config folder (required for some plugins, such as Condel
    mkdir -p /data/.vep/Plugins/config
    downloads+=("nxf_s3_retry nxf_s3_download s3://vep-databases/Plugins/config/ /data/.vep/Plugins/config/")

    # download annotation databases
    for ind in \${!annotation_clouds[@]}
    do
      cloud=\$(printf "\${annotation_clouds[\$ind]}" "\$chrom")
      localFile=\${annotation_locals[\$ind]}
      mkdir -p \$(dirname \${localFile})
      if [[ "\$cloud" == *gz ]]
      then
        mkdir -p \$(dirname \${localFile})
        if [[ "\$cloud" == gs* ]]
        then
          echo "this should be gnomad"
          downloads+=("/root/google-cloud-sdk/bin/gsutil cp \${cloud} \$localFile")
          downloads+=("/root/google-cloud-sdk/bin/gsutil cp \${cloud}.tbi \${localFile}.tbi")
        else
          downloads+=("nxf_s3_retry nxf_s3_download \$cloud \$localFile")
          downloads+=("nxf_s3_retry nxf_s3_download \${cloud}.tbi \${localFile}.tbi")
        fi
      else
        downloads+=("nxf_s3_retry nxf_s3_download \$cloud \$localFile")
      fi
    done

    # download plugin script
    mkdir -p /data/.vep/Plugins
    aws_profile="${params.vep_s3_profile}"
    downloads+=("nxf_s3_retry nxf_s3_download ${params.vep_plugins_path} /data/.vep/Plugins")

    # download cache
    # download relevant files in vep
    vep_files=\$(/home/ec2-user/miniconda/bin/aws s3 ls ${params.vep_s3_profile} ${params.vep_database_path}/${params.vep_cache}/\${chrom}/ | awk '{print \$4'})
    for vep_file in \${vep_files[@]}
    do
      is_overlap=\$(overlap \$start \$end \${vep_file})
      if [ \${is_overlap} -eq 1 ]; then
          source="${params.vep_database_path}/${params.vep_cache}/\${chrom}/\$vep_file"
          target="/data/.vep/${params.vep_cache}/\${chrom}/\$vep_file"
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
    mkdir -p "/data/.vep/${params.vep_cache}"
    for include in \${inclusion[@]}
    do
      downloads+=("nxf_s3_retry nxf_s3_download ${params.vep_database_path}/\${include} /data/.vep/\${include}")
    done

    # download plugin databases
    for ind in \${!plugin_clouds[@]}
    do
      cloud=\${plugin_clouds[\$ind]}
      localFile=\${plugin_locals[\$ind]}
      if [[ "\$cloud" == *gz ]]
      then
        mkdir -p \$(dirname \${localFile})
        if [[ "\$cloud" == gs* ]]
        then
          downloads+=("/root/google-cloud-sdk/bin/gsutil cp \${cloud} \$localFile")
          downloads+=("/root/google-cloud-sdk/bin/gsutil cp \${cloud}.tbi \${localFile}.tbi")
        else
          downloads+=("nxf_s3_retry nxf_s3_download \$cloud \$localFile")
          downloads+=("nxf_s3_retry nxf_s3_download \${cloud}.tbi \${localFile}.tbi")
        fi
      elif [ ! -z "\$cloud" ]
      then
        mkdir -p \$(dirname \${localFile})
        downloads+=("nxf_s3_retry nxf_s3_download \$cloud \$localFile")
      fi
    done

    # shuffle downloads since tabix will struggle
    #declare -a outarray
    #perm "\${downloads[@]}"

    nxf_parallel "\${downloads[@]}"

    # a plugin uses Math::CDF that is not in place
    cpanm Math::CDF

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
    plugins="\${plugins} --plugin CADD,${input_cadd}"

    annotations=
    for ind in \${!annotation_clouds[@]}
    do
      localfile=\${annotation_locals[\$ind]}
      cmd="zcat \${localfile} | grep -v '#' -c"
      #any_entry=\$(eval "\$cmd")
      any_entry=\$(if_has_record \$localfile)
      [[ \${any_entry} == 1 ]] && annotations="\${annotations} --custom \${localfile},\${annotation_annotations[\$ind]}"
    done
    annotations="\${annotations} --custom ${input_vcf},input,vcf,exact,0,${params.SELF_INFO_FIELDS}"
    ls -l /data/.vep/dbs/gnomad/

    vep -i ${input_vcf} ${params.vep_flags} --fork ${params.vep_threads} --dir /data/.vep --dir_cache /data/.vep --fasta ${human_ref} -a ${params.build} \$annotations \$plugins -o ${targetId}.vep.json
    

    """

}
/* somehow not working
process 'merge' {
  memory '7 G'
  cpus 1
  //label 'small_batch'
  //container params.picard_docker
  echo true
  publishDir ".", mode: 'copy'
  input:
    path(input_vcfs) from Cadd_ch.collect()
    path(input_indices) from Cadd_index_ch.collect()
  output:
    path("${cohortId}.cadd.tsv.gz*") into merge_ch

  
  """
  sample_list=(${input_vcfs})
  samples=\$(printf " I=%s " "\${sample_list[@]}")
  ls -alh

  java -jar /home/pontikos/bin/picard.jar MergeVcfs O=${cohortId}.cadd.tsv.gz \$samples
  touch ${cohortId}.cadd.tsv.gz
  """

}
*/
