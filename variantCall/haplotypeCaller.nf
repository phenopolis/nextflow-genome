/*
use haplotypeCaller to produce gVCFs
The input might come from another profile, so need to provide a list
NXF_VER=20.01.0 nextflow run variantCall.nf -w s3://phenopolis-nextflow/nextflow-work-test -with-report report.html
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

// calculate sample size
// since bam will come with bam.bai, will use `fromFilePairs`
// Even if Bam_ch will be used once by condition, Nextflow still complains for multiple use of the channel.

Calling_interval_list_ch = Channel.value(file("${calling_interval_list}"))

process 'getSamples' {
  output:
    path('file_list.txt') into Bam_file_list_ch

  """
  aws s3 ls ${params.input_bam_profile} ${params.input_bam_path}/ | awk '/${params.input_bam_suffix}\$/ {print \$4}' > file_list.txt
  """
}

Bam_ch = Bam_file_list_ch.splitText().map { it -> it.trim() }.map { it -> tuple(it - ~/${params.input_bam_suffix}$/, it) }

process 'HaplotypeCallerGvcf' {
  tag "$sampleId"
  label 'small_batch'
  memory '12 GB'
  container params.gatk_docker
  publishDir 's3://phenopolis-nextflow/gvcf', mode: 'copy'
  input:
    tuple val(sampleId), val(fileName) from Bam_ch
    file '*' from Calling_interval_list_ch
  output:
    tuple sampleId, path("${sampleId}.g.vcf.gz*") into HaplotypeCallerGvcf_VQSR_ch
  
  script:
    """
    # copy human_ref
    /home/ec2-user/miniconda/aws s3 cp ${params.aws_ref_profile} ${params.human_ref_path} . --recursive --exclude "*" --include "${params.human_ref_base}*"

    # copy sample
    /home/ec2-user/miniconda/aws s3 cp ${params.input_bam_profile} ${params.input_bam_path}/${fileName}.bam .
    /home/ec2-user/miniconda/aws s3 cp ${params.input_bam_profile} ${params.input_bam_path}/${fileName}.bam.bai .

    # gatk
    gatk --java-options \"${params.gatk_options} -Xmx6G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10\" \
      HaplotypeCaller \
        -R ${human_ref} \
        -I ${sampleId}.bam \
        -L ${exome_interval_list} \
        -O ${sampleId}.g.vcf.gz \
        -ERC GVCF
    """
}
