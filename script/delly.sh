#!/bin/bash

set -e
set -o pipefail

echo "********************"
echo "******Run delly*****"
echo "********************"


while getopts ":n:c:t:i:r:b:1:2:z:x:s:" opt; do
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
    r)
      refdata=$(realpath "$OPTARG")
      echo "reference: $refdata"
      ;;
    b)
      refbuild="$OPTARG"
      echo "build: $refbuild"
      ;;
    1)
      prefix1="$OPTARG"
      echo "prefix1: $prefix1"
      ;;
    2)
      prefix2="$OPTARG"
      echo "prefix2: $prefix2"
      ;;
    z)
      cnvsize="$OPTARG"
      echo "min CNV size:$cnvsize"
      ;;
    x)
      sdrd="$OPTARG"
      echo "min SD read-depth shift:$sdrd"
      ;;
    s)
      cnoffset="$OPTARG"
      echo "min CN offset:$cnoffset"
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
if [[ -z "$normal" || -z "$tumor" || -z "$refdata" || -z "$refbuild" || -z "$prefix1" || -z "$prefix2" || -z "$cnvsize" || -z "$sdrd" || -z "$cnoffset" ]]; then
  echo "One or more required arguments are missing." >&2
  exit 1
fi

workdir="${WORK_DIR}/${id}/delly"



if [[ -d "$workdir" ]]; then
  rm -rf $workdir
fi

mkdir $workdir
cd "$workdir" || exit

echo "${prefix1}"$'\t'"tumor" > samples.tsv
echo "${prefix2}"$'\t'"control" >> samples.tsv

if [ "$refbuild" == "hg19" ]; then
  map="${DATABASE_DIR}/database/delly_map/Homo_sapiens.GRCh37.dna.primary_assembly.fa.r101.s501.blacklist.gz"
fi

if [ "$refbuild" == "hg38" ]; then
  map="${DATABASE_DIR}/database/delly_map/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz"
fi

delly call -x ${TOOL_DIR}/delly/excludeTemplates/human.${refbuild}.excl.tsv  -o delly.sv.somatic.bcf -g $refdata $tumor $normal || { echo "Error during delly call sv bcf"; exit 1; }
delly filter -f somatic -o delly.sv.somatic.pre.bcf -s samples.tsv delly.sv.somatic.bcf || { echo "Error during delly filter"; exit 1; }
bcftools convert delly.sv.somatic.pre.bcf -o delly.sv.somatic.pre.vcf || { echo "Error during bcftools convert"; exit 1; }


delly cnv -u -z $cnvsize -x $sdrd -t $cnoffset -o delly.tumor.cnv.bcf -c delly.tumor.cov.gz -g $refdata -m $map $tumor || { echo "Error during delly call cnv bcf"; exit 1; }
delly cnv -u -x $sdrd -t $cnoffset -v delly.tumor.cnv.bcf -o delly.control.cnv.bcf -g $refdata -m $map $normal || { echo "Error during delly call cnv bcf"; exit 1; }

bcftools merge -m id -O b -o delly.tumor_control.cnv.bcf delly.tumor.cnv.bcf delly.control.cnv.bcf || { echo "Error during bcftools merge"; exit 1; }
bcftools index delly.tumor_control.cnv.bcf || { echo "Error during bcftools index"; exit 1; }
delly classify -p -f somatic -o delly.somatic.cnv.bcf -s samples.tsv delly.tumor_control.cnv.bcf || { echo "Error during delly cnv classify"; exit 1; }
bcftools convert delly.somatic.cnv.bcf -o delly.somatic.cnv.vcf || { echo "Error during bcftools convert"; exit 1; }
bcftools convert delly.tumor_control.cnv.bcf -o delly.tumor_control.cnv.vcf || { echo "Error during bcftools convert"; exit 1; }
bcftools query -s $prefix1 -f "%CHROM\t%POS\t%INFO/END\t%ID\t[%RDCN]\n" delly.somatic.cnv.vcf > segmentation.bed

echo "********************"
echo "*****Delly done*****"
echo "********************"
