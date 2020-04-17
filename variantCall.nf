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
    calling_interval_list = $params.exome_interval_list_hg38
  } else if (params.mode == 'wgs') {
    calling_interval_list = $params.wgs_interval_list_hg38
  }
  CNN_resources = params.CNN_resources_hg38
} else if (HG19.contains(params.build)) {
  // hg19
  human_ref_base = params.human_ref_base_hg19
  gatk_bundle = params.gatk_bundle_hg19

  // interval list for wes/wgs
  if (params.mode == 'wes'){
    calling_interval_list = $params.exome_interval_list_hg19
  } else if (params.mode == 'wgs') {
    calling_interval_list = $params.wgs_interval_list_hg19
  }
  CNN_resources = params.CNN_resources_hg19
} else {
  exit 1, "cannot understand build. Please provide either --build hg19, or --build hg38 (default is hg19)"
}
human_ref = "${human_ref_base}.fasta"

// calculate sample size
// since bam will come with bam.bai, will use `fromFilePairs`
sampleSize = Channel.fromFilePairs(params.input_bam).flatten().count().val / 3
Bam_ch = Channel.fromFilePairs(params.input_bam)
Human_ref_withdict_ch = Channel.value(file("${gatk_bundle}/${human_ref_base}*"))
Exome_interval_list_ch = Channel.value(file(exome_interval_list))

if (sampleSize < 30 && params.mode == 'wes')
  do_what = 'CNN'
else
  do_what = 'VQSR'


/********************************************************
 * CNN *
 * When sample size is small
 *******************************************************/
/*
 * part CNN_1 HaplotypeCallerGvcf for < 30 exomes
 */

process 'CNN_1_HaplotypeCallerGvcf' {
  cpus 3
  memory '16 G'
  container params.gatk_docker
  input:
    tuple val(sampleId), path(input_bam) from Bam_ch
    file '*' from Human_ref_withdict_ch
    file '*' from Exome_interval_list_ch
  output:
    tuple sampleId, path("hc.vcf.gz*"), path("bamout.bam*") into HaplotypeCallerGvcf_CNN_ch
  
  when:
    do_what == 'CNN'
  """
  gatk --java-options  \"${params.gatk_options} -Xmx7G\" \
    HaplotypeCaller \
      -R ${human_ref} \
      -I ${sampleId}.bam \
      -O hc.vcf.gz \
      -L ${calling_interval_list} \
      -bamout bamout.bam \
      --dont-trim-active-regions -stand-call-conf 0 -A Coverage -A ChromosomeCounts -A BaseQuality -A FragmentLength -A MappingQuality -A ReadPosition
  """
}

/*
 * part CNN_2 CNNScoreVariants for < 30 exomes
 */
process 'CNN_2_CNNScoreVariants' {
  cpus 3
  memory '8 G'
  container params.gatk_docker
  input:
    tuple val(sampleId), path(input_vcf), path(input_bam) from HaplotypeCallerGvcf_CNN_ch
    file '*' from Human_ref_withdict_ch
 output:
    tuple sampleId, path("CNN.vcf.gz*") into CNN_ch

  when:
    do_what == 'CNN'
  """
  gatk --java-options \"${params.gatk_options} -Xmx3G\" \
    CNNScoreVariants \
        -I bamout.bam \
        -R ${human_ref} \
        -V hc.vcf.gz \
        -O CNN.vcf.gz \
        -L ${calling_interval_list} \
        --tensor-type read_tensor \
        --inference-batch-size 8 \
        --transfer-batch-size 32 
  """
}

/*
 * part CNN_3 filter CNNScoreVariants for < 30 exomes
 */
GATK_ch = Channel.value(file("${gatk_bundle}/*"))
process 'CNN_3_filterCNNVariants' {
  cpus 3
  memory '20 G'
  container params.gatk_docker
  input:
    tuple val(sampleId), path(input_vcf) from CNN_ch
    file '*' from GATK_ch
  output:
    tuple sampleId, path("${sampleId}.CNN.filtered.vcf.gz*") into CNN_filtered_ch
  
  when:
    do_what == 'CNN'
  
  """
  # make resource list
  resource_list=(${params.CNN_resources_hg19})
  resource_list=\$(printf " --known-sites %s" "\${resource_list[@]}")
  resource_list=\${resource_list:1}
  gatk --java-options \"${params.gatk_options} -Xmx15G\" \
    FilterVariantTranches \
        -V CNN.vcf.gz \
        --output ${sampleId}.CNN.filtered.vcf.gz \
        \$resource_list \
        -info-key CNN_2D \
        --snp-tranche ${params.CNN_filter_snp_tranche} \
        --indel-tranche ${params.CNN_filter_indel_tranche}
  """
}

/*
 * part CNN_4 merge. Modify when to activate it. Need to see if filter is in genotype field before proceed.
 */
process 'CNN_4_merge_CNNVariants' {
  memory '20 G'
  cpus 3
  container params.gatk_docker
  publishDir "${params.outdir}/CNN_vcf", mode: 'copy'
  input:
    path(input_vcfs) from CNN_filtered_ch.collect()
  output:
    path("${params.cohort_name}.CNN.vcf.gz*") into CNN_merge_ch

  when:
    do_what == 'CNN'
  
  """
  bcftools merge -m all -Oz -o ${params.cohort_name}.CNN.vcf.gz ${input_vcfs}
  """
}

/****************************
 * VQSR
 * when sample size is big, or doing WGS
 ***************************/
/*
 *  part 2A HaplotypeCallerGvcf 
 */
process '2A_HaplotypeCallerGvcf' {
  memory '12 GB'
  container params.gatk_docker
  input:
    tuple val(sampleId), path(input_bam) from Bam_ch
    file '*' from Human_ref_withdict_ch
    file '*' from Exome_interval_list_ch
  output:
    tuple sampleId, path("${sampleId}.g.vcf.gz*") into HaplotypeCallerGvcf_VQSR_ch
  
  when:
    do_what == 'VQSR'

    """
    gatk --java-options \"${params.gatk_options} -Xmx6G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10\" \
      HaplotypeCaller \
        -R ${human_ref} \
        -I ${sampleId}.bam \
        -L ${exome_interval_list} \
        -O ${sampleId}.g.vcf.gz \
        -ERC GVCF
    """
}



/*
 * part 2B importGVCFs to genomicsDB
 if already exist, add to it
 */
process '2B_importGVCFs' {
  memory '8 G'
  container params.gatk_docker
  input:
    path(gvcfs) from haplotypeCallerGvcf_VQSR_ch.collect()

  output:
    val(true) into genomicsDB_ch

  when:
    do_what == 'VQSR'

  """
  # if output exist, add. else, create
  # !!!! -ip (padding) only applies to exome. remove it for genome
  ip=
  [[ \"${params.mode}\" == "wes" ]] && ip=\"-ip 500\"
  # make a sample string
  sample_list=(${gvcfs})
  samples=\$(printf " -V %s" "\${sample_list[@]}")
  samples=\${samples:1}

  if [ ! -d \"${params.genomicsDB}\" ] 
  then
    echo \"Create genomicsDB from scratch.\"
    gatk --java-options \"${params.gatk_options} -Xmx4g -Xms4g\" \
      GenomicsDBImport \
        --genomicsdb-workspace-path gendb://${params.genomicsDB} \
        --batch-size 5 \
        -L ${exome_interval_list} \
        -V \$samples \
        --reader-threads 5 \
        \$ip
  else
    echo \"genomicsDB already exists. Backup first, then add new gVCFs to it. https://gatk.broadinstitute.org/hc/en-us/articles/360035891051-GenomicsDB\"
    rm -rf ${params.results}/genomicsDB_backup
    cp -r ${params.results}/genomicsDB ${params.results}/genomicsDB_backup
    ${GATK} --java-options \"${params.gatk_options} -Xmx4g -Xms4g\"
      GenomicsDBImport \
        --genomicsdb-update-workspace-path gendb://${params.results}/genomicsDB \
        --batch-size 5 \
        -L ${params.interval_list} \
        -V \$samples \
        --reader-threads 5 \
        \$ip
  fi
  """
}

/*
 * part 2C GenotypeGVCFs
 */
process '2C_genotypeGVCFs' {
  input:
    val(flag) from genomicsDB_ch
  
  output:
    path("genotype.vcf.gz") into genotype_ch
  
  when:
    do_what == 'VQSR'
  
  """
  source /home/jingyu/.bashrc
  ${GATK} --java-options \"${params.gatk_options} -Xmx5g -Xms5g\"
    GenotypeGVCFs \
     -R ${params.genome} \
     -O genotype.vcf.gz \
     -D ${params.gatk_bundle}/dbsnp_138.b37.vcf.gz \
     -G StandardAnnotation \
     --only-output-calls-starting-in-intervals \
     --use-new-qual-calculator \
     -V gendb://${params.results}/genomicsDB \
     -L ${params.interval_list}
  """
}

/*
 * part 2D HardFilterAndMakeSitesOnlyVcf
 * Only does hard filter if do_excessHet is true
 */
process '2D_HardFilterAndMakeSitesOnlyVcf' {
  input:
    path(input_vcf) from genotype_ch
  
  output:
    path("siteOnly.vcf.gz") into SiteOnly_ch
  
  when:
    do_what == 'VQSR'

  script:
  if (params.do_excessHet)

  """
  source /home/jingyu/.bashrc
  echo \"Pipeline_logging: will do ExcessHet hard filtering\"
  ${GATK} --java-options \"${params.gatk_options} -Xmx3g -Xms3g\"
    VariantFiltration \
      --filter-expression "ExcessHet > ${params.excess_het_threshold}" \
      --filter-name ExcessHet \
      -O hardFiltered.vcf.gz \
      -V ${input_vcf} && \
  ${GATK} --java-options \"-Xmx3g -Xms3g\"
    MakeSitesOnlyVcf \
      --INPUT hardFiltered.vcf.gz \
      --OUTPUT siteOnly.vcf.gz
  """

  else

  """
  source /home/jingyu/.bashrc
  echo \"Pipeline_logging: will do not ExcessHet hard filtering\"

  ${GATK} --java-options \"${params.gatk_options} -Xmx3g -Xms3g\"
    MakeSitesOnlyVcf \
      --INPUT ${input_vcf} \
      --OUTPUT siteOnly.vcf.gz
  """
}

/*
 * part 2E IndelsVariantRecalibrator 
 */
process '2E_IndelsVariantRecalibrator' {
  memory params.mode=="wgs" ? 26.GB : 8.GB
  // clusterOptions completely override the one set in config, so need to restate all variables
  clusterOptions "-l h_rt=240:0:0 -S /bin/bash -l tmem=${memory.toString().replaceAll(/[\sB]/,'')},h_vmem=${memory.toString().replaceAll(/[\sB]/,'')}"
  input:
    path(input_vcf) from SiteOnly_ch

  output:
    tuple path("indel.recal"), path("indel.tranches") into Indel_recal_ch
  
  when:
    do_what == 'VQSR'
  
  script:
    def mem = params.mode == "wgs" ? "24g" : "6g"
    """
    source /home/jingyu/.bashrc
    ${GATK} --java-options \"${params.gatk_options} -Xmx${mem} -Xms${mem}\" \
    VariantRecalibrator \
      -V ${input_vcf} \
      -O indel.recal \
      --tranches-file indel.tranches \
      --trust-all-polymorphic \
      -tranche ${params.indel_recalibration_tranche_values} \
      -an ${indel_recalibration_annotation_values} \
      -mode INDEL \
      --max-gaussians 4 \
      --resource:mills,known=false,training=true,truth=true,prior=12 ${params.gatk_bundle}/Mills_and_1000G_gold_standard.indels.b37.vcf \
      --resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${params.gatk_bundle}/Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2 ${params.gatk_bundle}/dbsnp_138.b37.vcf.gz
    """
}

/*
 * part 2F SNPVariantRecalibrator create model when sample number >10000
 */
process '2F_SNPVariantRecalibrator_createModel' {
  memory 104.GB
  // clusterOptions completely override the one set in config, so need to restate all variables
  clusterOptions "-l h_rt=240:0:0 -S /bin/bash -l tmem=${memory.toString().replaceAll(/[\sB]/,'')},h_vmem=${memory.toString().replaceAll(/[\sB]/,'')}"
  input:
    path(input_vcf) from SiteOnly_ch
  
  output:
    path("snp.model.report") into SNP_model_ch
  when:
    do_what == 'VQSR' && sampleSize > 10000
  
  """
  source /home/jingyu/.bashrc
  ${GATK} --java-options \"${params.gatk_options} -Xmx100g -Xms100g\" \
    VariantRecalibrator \
      -V ${input_vcf} \
      -O  snp.recal \
      --tranches-file snp.tranches \
      --trust-all-polymorphic \
      -tranche ${params.SNP_recalibration_tranche_values} \
      -an ${SNP_recalibration_annotation_values} \
      -mode SNP \
      --sample-every-Nth-variant 10 \
      --output-model snp.model.report \
      --max-gaussians 6 \
      --resource:hapmap,known=false,training=true,truth=true,prior=15 ${params.gatk_bundle}/hapmap_3.3.b37.vcf \
      --resource:omni,known=false,training=true,truth=true,prior=12 ${params.gatk_bundle}/1000G_omni2.5.b37.vcf \
      --resource:1000G,known=false,training=true,truth=false,prior=10 ${params.gatk_bundle}/1000G_phase1.snps.high_confidence.b37.vcf.gz \
      --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${params.gatk_bundle}/dbsnp_138.b37.vcf.gz
  """

}

/*
 * part 2G SNPsVariantRecalibrator
 *  If dealing with many samples, might receive a model from previous step
 */
process '2G_SNPsVariantRecalibrator' {
  input:
    path(input_vcf) from SiteOnly_ch
    val model from SNP_model_ch.ifEmpty('EMPTY')
  output:
    tuple path("snp.recal"), path("snp.tranches") into SNP_recal_ch
  when:
    do_what == 'VQSR'
  
  script:
  def model_flag = model != 'EMPTY' ? "--input-model $model " : ''
  """
  source /home/jingyu/.bashrc
  ${GATK} --java-options \"${params.gatk_options} -Xmx3g -Xms3g\" \
    VariantRecalibrator \
      -V ${input_vcf} \
      -O snp.recal \
      --tranches-file snp.tranches \
      --trust-all-polymorphic \
      -tranche ${params.SNP_recalibration_tranche_values} \
      -an ${SNP_recalibration_annotation_values} \
      -mode SNP \
      ${model_flag} \
      --max-gaussians 6 \
      --resource:hapmap,known=false,training=true,truth=true,prior=15 ${params.gatk_bundle}/hapmap_3.3.b37.vcf \
      --resource:omni,known=false,training=true,truth=true,prior=12 ${params.gatk_bundle}/1000G_omni2.5.b37.vcf \
      --resource:1000G,known=false,training=true,truth=false,prior=10 ${params.gatk_bundle}/1000G_phase1.snps.high_confidence.b37.vcf.gz \
      --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${params.gatk_bundle}/dbsnp_138.b37.vcf.gz
  """
}

process '2H_ApplyRecalibration' {
  publishDir "${params.results}/vcf", mode: 'copy'
  input:
    path(input_vcf) from SiteOnly_ch
    tuple path(indel_recal), path(indel_tranch) from Indel_recal_ch
    tuple path(snp_recal), path(snp_tranch) from SNP_recal_ch
  output:
    path("${params.cohort_name}.vqsr.vcf") into VQSR_ch
  
  when:
    do_what == 'VQSR'

  """
  source /home/jingyu/.bashrc
  mkdir -p ${params.results}/vcf
  ${GATK} --java-options \"${params.gatk_options} -Xmx5g -Xms5g\" \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V ${input_vcf} \
      --recal-file ${indel_recal} \
      --tranches-file ${indel_tranch} \
      --truth-sensitivity-filter-level ${params.indel_filter_level} \
      --create-output-variant-index true \
      -mode INDEL

    ${GATK} --java-options \"${params.gatk_options} -Xmx5g -Xms5g\" \
      ApplyVQSR \
      -O ${params.cohort_name}.vqsr.vcf \
      -V tmp.indel.recalibrated.vcf \
      --recal-file ${snp_recal} \
      --tranches-file ${snp_tranch} \
      --truth-sensitivity-filter-level ${params.snp_filter_level} \
      --create-output-variant-index true \
      -mode SNP
  """
}