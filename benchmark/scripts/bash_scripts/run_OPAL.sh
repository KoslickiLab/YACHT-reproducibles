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

## go to the benchmark directory
cd $benchmark_dir
if [ ! -d $benchmark_dir/results/OPAL_results ]; then
    mkdir $benchmark_dir/results/OPAL_results
fi

# run OPAL on Rhizosphere challenge data
if [ ! -d $benchmark_dir/results/OPAL_results/rhizosphere_data ]; then
    mkdir $benchmark_dir/results/OPAL_results/rhizosphere_data;
fi
filename=`ls $benchmark_dir/CAMI_data/rhizosphere_data/ground_truth/*_sample_*.profile | awk -F '/' '{print $NF}' | head -1 | awk -F '_' '{print $1}'`
# run OPAL
parallel -j $cpu_num opal.py -g $benchmark_dir/CAMI_data/rhizosphere_data/ground_truth/${filename}_sample_{2}.profile -o $benchmark_dir/results/OPAL_results/rhizosphere_data/coverage{1}_sample{2} $benchmark_dir/results/YACHT_results/rhizosphere_data/${filename}_sample_{2}_cami_format_coverage{1}.profile ::: 1 0.5 0.1 0.05 0.01 ::: {0..20};
# Generate plots
python $benchmark_dir/scripts/python_scripts/generate_plots.py --cami_opal_res $benchmark_dir/CAMI_data/rhizosphere_data/opal_results/results.tsv --yacht_opal_dir $benchmark_dir/results/OPAL_results/rhizosphere_data --outdir $benchmark_dir/results/OPAL_results/rhizosphere_data --filename "rhizosphere" --fig_title "Rhizosphere Samples, Species Level"

# # run OPAL on Pathogen detection challenge data
# if [ ! -d $benchmark_dir/results/OPAL_results/pathogen_detection_data ]; then
#     mkdir $benchmark_dir/results/OPAL_results/pathogen_detection_data;
# fi
# filename=`ls $benchmark_dir/CAMI_data/pathogen_detection_data/ground_truth/*_sample_*.profile | awk -F '/' '{print $NF}' | head -1 | awk -F '_' '{print $1}'`
# # run OPAL
# parallel -j $cpu_num  opal.py -g $benchmark_dir/CAMI_data/pathogen_detection_data/ground_truth/${filename}_sample_{}.profile -o $benchmark_dir/results/OPAL_results/pathogen_detection_data/sample{} $benchmark_dir/results/YACHT_results/pathogen_detection_data/${filename}_sample_{}_cami_format.profile ::: {0..20};
# # Generate plots
# python $benchmark_dir/scripts/python_scripts/generate_plots.py --cami_opal_res $benchmark_dir/CAMI_data/pathogen_detection_data/opal_results/results.tsv --yacht_opal_dir $benchmark_dir/results/OPAL_results/rhizosphere_data --outfile $benchmark_dir/results/OPAL_results/pathogen_detection_data/pathogen_detection_comparison_plots.png


# run OPAL on Marine challenge data
if [ ! -d $benchmark_dir/results/OPAL_results/marine_data ]; then
    mkdir $benchmark_dir/results/OPAL_results/marine_data;
fi
filename=`ls $benchmark_dir/CAMI_data/marine_data/ground_truth/*_sample_*.profile | awk -F '/' '{print $NF}' | head -1 | awk -F '_' '{print $1}'`
# run OPAL
parallel -j $cpu_num opal.py -g $benchmark_dir/CAMI_data/marine_data/ground_truth/${filename}_sample_{2}.profile -o $benchmark_dir/results/OPAL_results/marine_data/coverage{1}_sample{2} $benchmark_dir/results/YACHT_results/marine_data/${filename}_sample_{2}_cami_format_coverage{1}.profile ::: 1 0.5 0.1 0.05 0.01 ::: {0..9};
# Generate plots
python $benchmark_dir/scripts/python_scripts/generate_plots.py --cami_opal_res $benchmark_dir/CAMI_data/marine_data/opal_results/results.tsv --yacht_opal_dir $benchmark_dir/results/OPAL_results/marine_data --outdir $benchmark_dir/results/OPAL_results/marine_data --filename "marine" --fig_title "Marine Samples, Species Level"

# run OPAL on Strain Madness challenge data
if [ ! -d $benchmark_dir/results/OPAL_results/strain_madness_data ]; then
    mkdir $benchmark_dir/results/OPAL_results/strain_madness_data;
fi
filename=`ls $benchmark_dir/CAMI_data/strain_madness_data/ground_truth/*_sample_*.profile | awk -F '/' '{print $NF}' | head -1 | awk -F '_' '{print $1}'`
# run OPAL
parallel -j $cpu_num opal.py -g $benchmark_dir/CAMI_data/strain_madness_data/ground_truth/${filename}_sample_{2}.profile -o $benchmark_dir/results/OPAL_results/strain_madness_data/coverage{1}_sample{2} $benchmark_dir/results/YACHT_results/strain_madness_data/${filename}_sample_{2}_cami_format_coverage{1}.profile ::: 1 0.5 0.1 0.05 0.01 ::: {0..99};
# Generate plots
python $benchmark_dir/scripts/python_scripts/generate_plots.py --cami_opal_res $benchmark_dir/CAMI_data/strain_madness_data/opal_results/results.tsv --yacht_opal_dir $benchmark_dir/results/OPAL_results/strain_madness_data --outdir $benchmark_dir/results/OPAL_results/strain_madness_data --filename "strain_madness" --fig_title "Strain Madness Samples, Species Level"

