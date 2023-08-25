#! /bin/bash
# run OPAL tool for YACHT algorithm 
# Usage:    

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: run_OPAL.sh <opal_repo_loc> <benchmark_dir> <cpu_num>"
    exit 1
fi
opal_repo_loc=$1
benchmark_dir=$2
cpu_num=$3


# # Split the ground truth profile into individual files
# cd $benchmark_dir/CAMI_data/rhizosphere_data/ground_truth
# awk '/@SampleID:/{filename=substr($0, length("@SampleID:") + 1) ".txt"} {print >filename}' gs_rhizosphere.profile
# # cd $benchmark_dir/CAMI_data/pathogen_detection_data/ground_truth
# # awk '/@SampleID:/{filename=substr($0, length("@SampleID:") + 1) ".txt"} {print >filename}' gs_pathogen.profile
# cd $benchmark_dir/CAMI_data/marine_data/ground_truth
# awk '/@SampleID:/{filename=substr($0, length("@SampleID:") + 1) ".txt"} {print >filename}' gs_marine_short.profile
# cd $benchmark_dir/CAMI_data/strain_madness_data/ground_truth
# awk '/@SampleID:/{filename=substr($0, length("@SampleID:") + 1) ".txt"} {print >filename}' gs_strain_madness_short_long.profile


## go to the benchmark directory
cd $benchmark_dir
if [ ! -d $benchmark_dir/results/OPAL_results ]; then
    mkdir $benchmark_dir/results/OPAL_results
fi

# run OPAL on Rhizosphere challenge data
if [ ! -d $benchmark_dir/results/OPAL_results/rhizosphere_data ]; then
    mkdir $benchmark_dir/results/OPAL_results/rhizosphere_data;
fi
cd $benchmark_dir/results/OPAL_results/rhizosphere_data;
filename=`ls $benchmark_dir/CAMI_data/rhizosphere_data/ground_truth/*.txt | awk -F '/' '{print $NF}' | head -1 | awk -F '_' '{print $1}'`
parallel -j $cpu_num python $opal_repo_loc/opal.py -g $benchmark_dir/CAMI_data/rhizosphere_data/ground_truth/${filename}_sample_{}.txt -o $benchmark_dir/results/OPAL_results/rhizosphere_data $benchmark_dir/results/YACHT_results/rhizosphere_data/${filename}_sample_{}_cami_format.profile ::: {0..20};