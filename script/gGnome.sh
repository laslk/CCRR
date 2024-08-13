#!/bin/bash

set -e
set -o pipefail

echo "******************"
echo "*****gGnome*******"
echo "*****************"
while getopts ":p:f:j:s:c:o:u:l:" opt; do
  case $opt in
    p)
      threads="$OPTARG"
      echo "threads: $threads"
      ;;
    j)
      vmheap="$OPTARG"
      echo "vm heap:$vmheap G"
      ;;
    f)
      config="$OPTARG"
      echo "config: $config"
      ;;
    s)
      sv="$OPTARG"
      echo "svfile: $sv"
      ;;
    c)
      cn="$OPTARG"
      echo "cnfile: $cn"
      ;;
    o)
      wd="$OPTARG"
      echo "wd: $wd"
      ;;
    u)
      cellularity="$OPTARG"
      echo "ploidy: $cellularity"
      ;;
    l)
      ploidy="$OPTARG"
      echo "ploidy: $ploidy"
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


export GUROBI_HOME=${TOOL_DIR}/jabba/gurobi1003/linux64
export LD_LIBRARY_PATH=${TOOL_DIR}/jabba/gurobi1003/linux64/lib
export PATH=/opt/conda/envs/main/bin:/opt/conda/condabin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:${TOOL_DIR}/jabba/gurobi1003/linux64/bin:/opt/conda/envs/main/lib/R/library/JaBbA/extdata

if [[ -d $wd ]]; then
  rm -rf $wd
fi
mkdir $wd
cd $wd

conda run -n main --no-capture-output Rscript ${SCRIPT_DIR}/format_conv/sv2gGnome.R $sv ${wd}/sv.bedpe || { echo "sv input Error"; exit 1; }
if [[ "$cellularity" != "null" && "$ploidy" != "null" ]]; then
  conda run -n main --no-capture-output jba ${wd}/sv.bedpe $cn --gurobi TRUE --verbose --outdir=${wd} --cores=$threads --mem=$vmheap --ploidy=$ploidy --purity=$cellularity || { echo "Error during jabba"; exit 1; }
else
  conda run -n main --no-capture-output jba ${wd}/sv.bedpe $cn --gurobi TRUE --verbose --outdir=${wd} --cores=$threads --mem=$vmheap || { echo "Error during jabba"; exit 1; }
fi


if [[ -d ${wd}/../gGnome ]]; then
  rm -rf ${wd}/../gGnome
fi
mkdir ${wd}/../gGnome
conda run -n main --no-capture-output Rscript ${SCRIPT_DIR}/gGnome.R ${wd} || { echo "Error during gGnome"; exit 1; }