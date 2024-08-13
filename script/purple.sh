#!/bin/bash

set -e
set -o pipefail

gridss_database="${DATABASE_DIR}/database/gridss_database/"

echo "********************"
echo "*****run purple*****"
echo "********************"


while getopts ":n:t:r:i:p:b:j:1:2:g:" opt; do
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
    g)
      using_gridss="$OPTARG"
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


if [[ -d "${WORK_DIR}/${id}/purple" ]]; then
  rm -rf ${WORK_DIR}/${id}/purple
fi
mkdir ${WORK_DIR}/${id}/purple
mkdir ${WORK_DIR}/${id}/purple/output
cd ${WORK_DIR}/${id}/purple

if [ "$using_gridss" == "True" ];then
    if [ "$ref" == "hg19" ];then
    java -Xmx${jvmheap}G -jar /opt/conda/envs/main/share/hmftools-purple-3.7.1-0/purple.jar -reference $name_n -tumor $name_t -amber ${WORK_DIR}/${id}/pre_purple/amber -cobalt ${WORK_DIR}/${id}/pre_purple/cobalt -gc_profile ${gridss_database}37/copy_number/GC_profile.1000bp.37.cnp -ref_genome $refdata -ref_genome_version 37 -ensembl_data_dir ${gridss_database}37/common/ensembl_data/ -somatic_vcf ${WORK_DIR}/${id}/pre_purple/sage/sage.vcf.gz -structural_vcf ${WORK_DIR}/${id}/gridss/gripss_output/${name_t}.gripss.filtered.vcf.gz -output_dir ${WORK_DIR}/${id}/purple/output/ -threads $threads || { rm -rf ${WORK_DIR}/${id}/purple;echo "Error during purple"; exit 1; }
    elif [ "$ref" == "hg38" ];then
    java -Xmx${jvmheap}G -jar /opt/conda/envs/main/share/hmftools-purple-3.7.1-0/purple.jar -reference $name_n -tumor $name_t -amber ${WORK_DIR}/${id}/pre_purple/amber -cobalt ${WORK_DIR}/${id}/pre_purple/cobalt -gc_profile ${gridss_database}38/copy_number/GC_profile.1000bp.38.cnp -ref_genome $refdata -ref_genome_version 38 -ensembl_data_dir ${gridss_database}38/common/ensembl_data/ -somatic_vcf ${WORK_DIR}/${id}/pre_purple/sage/sage.vcf.gz -structural_vcf ${WORK_DIR}/${id}/gridss/gripss_output/${name_t}.gripss.filtered.vcf.gz -output_dir ${WORK_DIR}/${id}/purple/output/ -threads $threads || { rm -rf ${WORK_DIR}/${id}/purple;echo "Error during purple"; exit 1; }
    fi
elif [ "$using_gridss" == "False" ];then
    if [ "$ref" == "hg19" ];then
    java -Xmx${jvmheap}G -jar /opt/conda/envs/main/share/hmftools-purple-3.7.1-0/purple.jar -reference $name_n -tumor $name_t -amber ${WORK_DIR}/${id}/pre_purple/amber -cobalt ${WORK_DIR}/${id}/pre_purple/cobalt -gc_profile ${gridss_database}37/copy_number/GC_profile.1000bp.37.cnp -ref_genome $refdata -ref_genome_version 37 -ensembl_data_dir ${gridss_database}37/common/ensembl_data/ -somatic_vcf ${WORK_DIR}/${id}/pre_purple/sage/sage.vcf.gz -output_dir ${WORK_DIR}/${id}/purple/output/ -threads $threads || { rm -rf ${WORK_DIR}/${id}/purple;echo "Error during purple"; exit 1; }
    elif [ "$ref" == "hg38" ];then
    java -Xmx${jvmheap}G -jar /opt/conda/envs/main/share/hmftools-purple-3.7.1-0/purple.jar -reference $name_n -tumor $name_t -amber ${WORK_DIR}/${id}/pre_purple/amber -cobalt ${WORK_DIR}/${id}/pre_purple/cobalt -gc_profile ${gridss_database}38/copy_number/GC_profile.1000bp.38.cnp -ref_genome $refdata -ref_genome_version 38 -ensembl_data_dir ${gridss_database}38/common/ensembl_data/ -somatic_vcf ${WORK_DIR}/${id}/pre_purple/sage/sage.vcf.gz -output_dir ${WORK_DIR}/${id}/purple/output/ -threads $threads || { rm -rf ${WORK_DIR}/${id}/purple;echo "Error during purple"; exit 1; }
    fi
fi

echo "*****************"
echo "***purple done***"
echo "*****************"
