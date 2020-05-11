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

Scattered_ch = Horizontal_split_ch
  .flatten()
  .map { file -> tuple(file.name.take(file.name.lastIndexOf('.vcf.gz')), file ) }
  .map { it -> tuple(it[0], file(it[1]), file("${params.input_beds}/${it[0]}.bed"))}

//cadd_files = file(params.cadd_file_list).readLines().join(' ')

process 'cadd' {
  tag "${targetId}"
  cpus 4
  label 'small_batch'
  memory '16 G'
  maxForks 1
  echo true

  container params.cadd_docker

  input:
    tuple val(targetId), path(input_vcf), path(input_bed) from Scattered_ch
  output:
    tuple targetId, path("${targetId}.cadd.tsv.gz") into Cadd_ch
    path("${targetId}.cadd.tsv.gz.tbi") into CaddIndex_ch

  """
  source s3.bash
  # copy CADD annotation locally using tabix for the files in the `files_tobe_tabix_download`,
  # and copy straight if not
  # copy to /cadd/data
  mkdir -p /cadd/data

  workdir=\$(pwd)
  # tabix
  chrom=\$(head -1 ${input_bed} | cut -f1)
  start=\$(head -1 ${input_bed} | cut -f2)
  start=\$(( \${start}-${params.input_padding} ))
  start=\$(( \${start} > 0 ? \${start} : 0 ))
  end=\$(tail -1 ${input_bed} | cut -f3)
  end=\$(( \${end}+${params.input_padding} ))
  tabix_files=(${params.cadd_files_tobe_tabix_download})
  aws_profile="${params.aws_cadd_profile}"
  # make exclusive list
  exclusion=\$(printf ' --exclude "%s*"' "\${tabix_files[@]}")
  exclusion=\${exclusion:1}
  echo \$exclusion
  nxf_s3_retry nxf_s3_download ${params.cadd_base_path} /cadd/data \$exclusion

  downloads=()
  for file in \${tabix_files[@]}
  do
    source=${params.cadd_base_path}/\${file}
    target=/cadd/data/\${file}
    target_dir=\$(dirname "\${target}")
    mkdir -p \$target_dir
    # some files have the same basename, which produces confusing error message at bgzip
    downloads+=("cd \$target_dir && HTS_S3_HOST=s3.eu-central-1.wasabisys.com AWS_PROFILE=wasabi tabix -h \$source \${chrom}:\${start}-\${end} | bgzip -c > \${target} && tabix -p vcf \${target} && cd \$workdir ")
  done

  nxf_parallel "\${downloads[@]}"

  cd \${workdir}

  # run cadd
  cadd ${params.cadd_flags} -o ${targetId}.cadd.tsv.gz ${input_vcf} && tabix -p vcf ${targetId}.cadd.tsv.gz
  """
}

