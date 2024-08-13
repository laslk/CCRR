#!/bin/bash

set -e
set -o pipefail


varscandata="${DATABASE_DIR}/database/varscan/"
varscanjar="/opt/conda/envs/main/share/varscan-2.4.6-0/VarScan.jar"

echo "********************"
echo "****Run varscan*****"
echo "********************"


while getopts ":n:t:j:r:p:b:i:" opt; do
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
      p="$OPTARG"
      echo "varscan somatic-p-value: $p"
      ;;
    b)
      ref="$OPTARG"
      echo "reference build: $ref"
      ;;
    j)
      jvmheap="$OPTARG"
      echo "jvm heap:$jvmheap G"
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


if [[ -d "${WORK_DIR}/${id}/varscan" ]]; then
  rm -rf  ${WORK_DIR}/${id}/varscan
fi

mkdir ${WORK_DIR}/${id}/varscan
cd "${WORK_DIR}/${id}/varscan"


samtools mpileup --no-BAQ -q 1 -f $refdata $normal $tumor > to_varscan.mpileup || { rm -rf ${WORK_DIR}/${id}/varscan;echo "Error during samtools mpileup"; exit 1; }

java -Xmx${jvmheap}G -jar $varscanjar somatic to_varscan.mpileup ${WORK_DIR}/${id}/varscan/varscan --mpileup 1 --output-vcf  --somatic-p-value $p || { rm -rf ${WORK_DIR}/${id}/varscan; echo "Error during varscan"; exit 1; }

java -Xmx${jvmheap}G -jar $varscanjar copynumber ${WORK_DIR}/${id}/varscan/to_varscan.mpileup varscan -mpileup 1 || { rm -rf ${WORK_DIR}/${id}/varscan; echo "Error during varscan copynumber"; exit 1; }

java -Xmx${jvmheap}G -jar $varscanjar copyCaller ${WORK_DIR}/${id}/varscan/varscan.copynumber --output-file varscan.copynumber.called || { rm -rf ${WORK_DIR}/${id}/varscan; echo "Error during varscan copyCaller"; exit 1; }

rm -f ${WORK_DIR}/${id}/varscan/to_varscan.mpileup
java -Xmx${jvmheap}G -jar $varscanjar processSomatic varscan.indel.vcf
java -Xmx${jvmheap}G -jar $varscanjar processSomatic varscan.snp.vcf

cat varscan.indel.Somatic.vcf varscan.snp.Somatic.vcf > varscan.somatic.vcf


echo "********************"
echo "***Varscan2 done****"
echo "********************"

