#!/bin/bash

set -e
set -o pipefail

echo "********************"
echo "****Run sequenza****"
echo "********************"


while getopts ":n:t:i:r:d:" opt; do
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
    d)
      threads="$OPTARG"
      echo "threads: $threads"
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

workdir="${WORK_DIR}/${id}/sequenza"

if [[ -d "$workdir" ]]; then
  rm -rf $workdir
fi
mkdir $workdir
cd "$workdir" || exit


refdir=$(dirname $refdata)
wig_file="${refdir}/$(basename $refdata).gc50Base.wig.gz"
if [ -e "${wig_file}" ]; then
  echo "File ${wig_file} already exists."
else
  sequenza-utils gc_wiggle -w 50 --fasta ${refdata} -o ${wig_file}
  echo "Generated file: ${wig_file}"
fi

sequenza-utils bam2seqz -n ${normal} -t ${tumor} --fasta ${refdata} -gc ${wig_file} -o seqz.gz
sequenza-utils seqz_binning --seqz ${workdir}/seqz.gz -w 50 > ${workdir}/small.seqz
Rscript ${SCRIPT_DIR}/sequenza.R $workdir $id $threads
rm -rf seqz.gz

echo "********************"
echo "***sequenza done****"
echo "********************"