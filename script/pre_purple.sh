#!/bin/bash

gridss_database="${DATABASE_DIR}/database/gridss_database/"

echo "********************"
echo "*prepare for purple*"
echo "********************"


while getopts ":n:t:r:p:b:j:1:2:i:" opt; do
  case $opt in
    n)
      normal=$(realpath "$OPTARG")
      echo "normal bam: $normal"
      ;;
    t)
      tumor=$(realpath "$OPTARG")
      echo "tumor bam: $tumor"
      ;;
    r)
      refdata=$(realpath "$OPTARG")
      echo "reference: $refdata"
      ;;
    i)
      id="$OPTARG"
      ;;
    p)
      threads="$OPTARG"
      echo "threads: $threads"
      ;;
    b)
      ref="$OPTARG"
      echo "reference build: $ref"
      ;;
    j)
      jvmheap="$OPTARG"
      echo "jvm heap:$jvmheap G"
      ;;
    1)
      name_n="$OPTARG"
      echo "normal name:$name_n"
      ;;
    2)
      name_t="$OPTARG"
      echo "tumor name:$name_t"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [[ -d "${WORK_DIR}/${id}/pre_purple" ]]; then
  rm -rf ${WORK_DIR}/${id}/pre_purple
fi
mkdir ${WORK_DIR}/${id}/pre_purple


if [[ ! -d "${WORK_DIR}/${id}/pre_purple/cobalt" ]]; then
  mkdir ${WORK_DIR}/${id}/pre_purple/cobalt
fi
cd ${WORK_DIR}/${id}/pre_purple/cobalt


if [ "$ref" == "hg19" ];then
  java -Xmx${jvmheap}G -jar /opt/conda/envs/main/share/hmftools-cobalt-1.13-1/cobalt.jar -reference $name_n -reference_bam $normal -tumor $name_t -tumor_bam $tumor -output_dir ${WORK_DIR}/${id}/pre_purple/cobalt -threads $threads -gc_profile ${gridss_database}37/copy_number/GC_profile.1000bp.37.cnp || { rm -rf ${WORK_DIR}/${id}/pre_purple/cobalt;echo "Error during cobalt"; exit 1; }
elif [ "$ref" == "hg38" ];then
  java -Xmx${jvmheap}G -jar /opt/conda/envs/main/share/hmftools-cobalt-1.13-1/cobalt.jar -reference $name_n -reference_bam $normal -tumor $name_t -tumor_bam $tumor -output_dir ${WORK_DIR}/${id}/pre_purple/cobalt -threads $threads -gc_profile ${gridss_database}38/copy_number/GC_profile.1000bp.38.cnp || { rm -rf ${WORK_DIR}/${id}/pre_purple/cobalt;echo "Error during cobalt"; exit 1; }
fi
echo "***COBALT done***"



if [[ ! -d "${WORK_DIR}/${id}/pre_purple/amber" ]]; then
  mkdir ${WORK_DIR}/${id}/pre_purple/amber
fi

cd ${WORK_DIR}/${id}/pre_purple/amber

if [ "$ref" == "hg19" ];then
  java -Xmx${jvmheap}G -cp /opt/conda/envs/main/share/hmftools-amber-3.9-1/amber.jar com.hartwig.hmftools.amber.AmberApplication -reference $name_n -reference_bam $normal -tumor $name_t -tumor_bam $tumor -ref_genome $refdata -output_dir ${WORK_DIR}/${id}/pre_purple/amber -threads $threads -ref_genome_version HG19 -loci ${gridss_database}37/copy_number/GermlineHetPon.37.vcf.gz || { rm -rf ${WORK_DIR}/${id}/pre_purple/amber;echo "Error during amber"; exit 1; }
elif [ "$ref" == "hg38" ];then
  java -Xmx${jvmheap}G -cp /opt/conda/envs/main/share/hmftools-amber-3.9-1/amber.jar com.hartwig.hmftools.amber.AmberApplication -reference $name_n -reference_bam $normal -tumor $name_t -tumor_bam $tumor -ref_genome $refdata -output_dir ${WORK_DIR}/${id}/pre_purple/amber -threads $threads -ref_genome_version V38 -loci ${gridss_database}38/copy_number/GermlineHetPon.38.vcf.gz || { rm -rf ${WORK_DIR}/${id}/pre_purple/amber;echo "Error during amber"; exit 1; }
fi
echo "***AMBER done***"


if [[ ! -d "${WORK_DIR}/${id}/pre_purple/sage" ]]; then
  mkdir ${WORK_DIR}/${id}/pre_purple/sage
fi
cd ${WORK_DIR}/${id}/pre_purple/sage


if [ "$ref" == "hg19" ];then
  java -Xms4G -Xmx${jvmheap}G -cp /opt/conda/envs/main/share/hmftools-sage-3.2.3-0/sage.jar com.hartwig.hmftools.sage.SageApplication -tumor $name_t -tumor_bam $tumor -ref_genome $refdata -ref_genome_version 37 -threads $threads -reference $name_n -reference_bam $normal -hotspots ${gridss_database}37/variants/KnownHotspots.somatic.37.vcf.gz  -panel_bed ${gridss_database}37/variants/ActionableCodingPanel.37.bed.gz -high_confidence_bed ${gridss_database}37/variants/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz -ensembl_data_dir ${gridss_database}37/common/ensembl_data/ -out sage.vcf.gz || { rm -rf ${WORK_DIR}/${id}/pre_purple/sage;echo "Error during sage"; exit 1; }
elif [ "$ref" == "hg38" ];then
  java -Xms4G -Xmx${jvmheap}G -cp /opt/conda/envs/main/share/hmftools-sage-3.2.3-0/sage.jar com.hartwig.hmftools.sage.SageApplication -tumor $name_t -tumor_bam $tumor -ref_genome $refdata -ref_genome_version 38 -threads $threads -reference $name_n -reference_bam $normal -hotspots ${gridss_database}38/variants/KnownHotspots.somatic.38.vcf.gz  -panel_bed ${gridss_database}38/variants/ActionableCodingPanel.38.bed.gz -high_confidence_bed ${gridss_database}38/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz -ensembl_data_dir ${gridss_database}38/common/ensembl_data/ -out sage.vcf.gz || { rm -rf ${WORK_DIR}/${id}/pre_purple/sage;echo "Error during sage"; exit 1; }
fi
echo "***SAGE done***"

