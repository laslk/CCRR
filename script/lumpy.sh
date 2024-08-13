#!/bin/bash

set -e
set -o pipefail


echo "********************"
echo "******Run Lumpy*****"
echo "********************"


while getopts ":n:t:i:" opt; do
  case $opt in
    n)
      normal=$(realpath "$OPTARG")
      echo "normal bam: $normal"
      ;;
    t)
      tumor=$(realpath "$OPTARG")
      echo "tumor bam: $tumor"
      ;;
    i)
      id="$OPTARG"
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


if [[ -d "${WORK_DIR}/${id}/lumpy" ]]; then
  rm -rf ${WORK_DIR}/${id}/lumpy
  mkdir ${WORK_DIR}/${id}/lumpy
elif [[ ! -d "${WORK_DIR}/${id}/lumpy" ]];then
  mkdir ${WORK_DIR}/${id}/lumpy
fi


mkdir ${WORK_DIR}/${id}/lumpy/pre_lumpy
cd ${WORK_DIR}/${id}/lumpy/pre_lumpy


samtools view -b -F 1294 $normal > normal.discordants.unsorted.bam || { echo "Error during samtools"; exit 1; }
samtools view -b -F 1294 $tumor > tumor.discordants.unsorted.bam || { echo "Error during samtools"; exit 1; }

samtools view -h normal.discordants.unsorted.bam | extractSplitReads_BwaMem -i stdin | samtools view -Sb - > normal.splitters.unsorted.bam || { echo "Error during samtools"; exit 1; }
samtools view -h tumor.discordants.unsorted.bam | extractSplitReads_BwaMem -i stdin | samtools view -Sb - > tumor.splitters.unsorted.bam || { echo "Error during samtools"; exit 1; }

samtools sort normal.discordants.unsorted.bam -o normal.discordants.bam || { echo "Error during samtools"; exit 1; }
samtools sort normal.splitters.unsorted.bam -o normal.splitters.bam || { echo "Error during samtools"; exit 1; }

samtools sort tumor.discordants.unsorted.bam -o tumor.discordants.bam || { echo "Error during samtools"; exit 1; }
samtools sort tumor.splitters.unsorted.bam -o tumor.splitters.bam || { echo "Error during samtools"; exit 1; }

rm -rf  *.unsorted.bam

cd ${WORK_DIR}/${id}/lumpy/

lumpyexpress -B $tumor,$normal -S pre_lumpy/tumor.splitters.bam,pre_lumpy/normal.splitters.bam -D pre_lumpy/tumor.discordants.bam,pre_lumpy/normal.discordants.bam -o lumpy.tumor_normal.vcf || { echo "Error during lumpy"; exit 1; }

svtyper -B $tumor,$normal -i lumpy.tumor_normal.vcf > lumpy.gt.vcf
rm -rf pre_lumpy

echo "********************"
echo "*****lumpy done*****"
echo "********************"

