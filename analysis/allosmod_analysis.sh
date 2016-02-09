#!/bin/bash

analysis_dir='/netapp/sali/peterc/cryptosite/src_multichain/analysis'
#if test -z $analysis_dir; then echo "set analysis_dir variable in allosmod_analysis.sh"; exit; fi

echo "get_e"
${analysis_dir}/bin/get_e.sh
echo "get_p"
${analysis_dir}/bin/get_p.sh
echo "get_v"
${analysis_dir}/bin/get_velocities.sh $analysis_dir
echo "get_qi"
${analysis_dir}/bin/get_qi.sh $analysis_dir

${analysis_dir}/bin/check_runs.sh >check_runs.out
cat check_runs.out

