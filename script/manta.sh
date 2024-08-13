#!/bin/bash

set -e
set -o pipefail



echo "********************"
echo "*****Run manta******"
echo "********************"


while getopts ":n:t:r:g:p:i:" opt; do
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
      echo "number of jobs: $threads"
      ;;
    g)
      memgb="$OPTARG"
      echo "available memory: ${memgb} g"
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


workdir="${WORK_DIR}/${id}/manta"

if [[ -d "$workdir" ]]; then
  rm -rf $workdir
fi

mkdir $workdir
cd "$workdir" || exit


configManta.py --normalBam $normal --tumorBam $tumor --referenceFasta $refdata --runDir $workdir || { rm -rf $workdir;echo "Error during manta"; exit 1; }
${workdir}/runWorkflow.py -j $threads -g $memgb || { rm -rf $workdir;echo "Error during manta"; exit 1; }


echo "********************"
echo "*****manta done*****"
echo "********************"
