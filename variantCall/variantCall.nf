/*
Variant call: from Bam to un-annotated vcf (with filter)
Note, it uses CNNScoreVariant
NXF_VER=20.01.0 nextflow run variantCall_CNN.nf --input_table input_bam.tsv -bucket-dir s3://phenopolis-nextflow/nextflow-work/variantCallCNN -with-report report.html
*/
// determines build
HG38 = ['hg38', 'GRCh38', 'grch38']
HG19 = ['hg19', 'GRCh37', 'b37', 'grch37']
human_ref_base = null
gatk_bundle = null
calling_interval_list = null
CNN_resources = null
if (HG38.contains(params.build)) {
  // hg38
  human_ref_base = params.human_ref_base_hg38
  gatk_bundle = params.gatk_bundle_hg38
  // interval list for wes/wgs
  if (params.mode == 'wes'){
    calling_interval_list = params.exome_interval_list_hg38
  } else if (params.mode == 'wgs') {
    calling_interval_list = params.wgs_interval_list_hg38
  }
  CNN_resources = params.CNN_resources_hg38
} else if (HG19.contains(params.build)) {
  // hg19
  human_ref_base = params.human_ref_base_hg19
  gatk_bundle = params.gatk_bundle_hg19

  // interval list for wes/wgs
  if (params.mode == 'wes'){
    calling_interval_list = params.exome_interval_list_hg19
  } else if (params.mode == 'wgs') {
    calling_interval_list = params.wgs_interval_list_hg19
  }
  CNN_resources = params.CNN_resources_hg19
} else {
  exit 1, "cannot understand build. Please provide either --build hg19, or --build hg38 (default is hg19)"
}
human_ref = "${human_ref_base}.fasta"
human_ref_all_bundle = "'${human_ref}' '${human_ref}.fai' '${human_ref}.pac' '${human_ref_base}.dict' '${human_ref}.amb' '${human_ref}.ann' '${human_ref}.bwt' '${human_ref}.sa'"

// Even if Bam_ch will be used once by condition, Nextflow still complains for multiple use of the channel.
Bam_ch = Channel.fromList(file(params.input_table).readLines())
interval_ch = Channel.value(file(calling_interval_list))


process 'CNN_HaplotypeCallerGvcf' {
  tag "$sampleId"
  cpus 4
  memory '16 G'
  label 'gatk'
  container params.gatk_docker
  input:
    val sampleId from Bam_ch
    file interval from interval_ch
  output:
    tuple sampleId, path("hc.vcf.gz*"), path("bamout.ba*") into HaplotypeCallerGvcf_CNN_ch
  
  """
  source s3.bash

  aws_profile="${params.s3_deposit_profile}"
  downloads=()
  input_bam=
  for filename in \$(/home/ec2-user/miniconda/bin/aws s3 ls ${params.s3_deposit_profile} ${params.s3_deposit}/${sampleId}/ | grep -e 'bam\$\\|bai\$' | grep -v 'disc_sorted.bam' | awk '{print \$NF}'); do
      cloudFile=${params.s3_deposit}/${sampleId}/\$filename
      downloads+=("nxf_s3_retry nxf_s3_download \$cloudFile \$filename")
      if [[ \$filename == *.bam ]]; then
        input_bam=\$filename
      fi
  done
  nxf_parallel "\${downloads[@]}"
  gatk --java-options  \"${params.gatk_options} -Xmx7G\" \
    HaplotypeCaller \
      -R /data/gatk/${human_ref} \
      -I \${input_bam} \
      -O hc.vcf.gz \
      -L ${interval} \
      -bamout bamout.bam \
      --dont-trim-active-regions -stand-call-conf 0 -A Coverage -A ChromosomeCounts -A BaseQuality -A FragmentLength -A MappingQuality -A ReadPosition
  """
}

process 'CNNScoreVariants' {
  tag "$sampleId"
  cpus 4
  memory '8 G'
  label 'gatk'
  container params.gatk_docker
  input:
    tuple val(sampleId), path(input_vcf), path(input_bam) from HaplotypeCallerGvcf_CNN_ch
    file interval from interval_ch
  output:
    tuple sampleId, path("CNN.vcf.gz*") into CNN_ch

  """
  samtools index bamout.bam
  gatk --java-options \"${params.gatk_options} -Xmx3G\" \
    CNNScoreVariants \
        -I bamout.bam \
        -R /data/gatk/${human_ref} \
        -V hc.vcf.gz \
        -O CNN.vcf.gz \
        -L ${interval} \
        --tensor-type read_tensor \
        --inference-batch-size 8 \
        --transfer-batch-size 32 
  """
}

process 'filterCNNVariants' {
  tag "$sampleId"
  cpus 1
  memory '20 G'
  publishDir "${params.publish_dir}/${sampleId}"
  label 'gatk'
  container params.gatk_docker
  input:
    tuple val(sampleId), path(input_vcf) from CNN_ch
  output:
    path("CNN.filtered.vcf.gz*") into Result_ch
  
  """
  # make resource list
  resource_list=(${params.CNN_resources_hg19})
  resource_list=\$(printf " -resource /data/gatk/%s" "\${resource_list[@]}")
  gatk --java-options \"${params.gatk_options} -Xmx15G\" \
    FilterVariantTranches \
        -V CNN.vcf.gz \
        --output CNN.filtered.vcf.gz \
        \$resource_list \
        -info-key CNN_2D \
        --snp-tranche ${params.CNN_filter_snp_tranche} \
        --indel-tranche ${params.CNN_filter_indel_tranche}
  """
}
