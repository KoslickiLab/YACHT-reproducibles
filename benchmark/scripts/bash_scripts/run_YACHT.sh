#! /bin/bash
# run benchmarking experiments based on CAMI2 data for YACHT algorithm
# Usage:    

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: run_YACHT.sh <yacht_repo_loc> <benchmark_dir> <cpu_num>"
    exit 1
fi
yacht_repo_loc=$1
benchmark_dir=$2
cpu_num=$3

## Run YACHT on CAMI2 data (e.g., Rhizosphere challenge, Clinical pathogen detection challenge, Challenge Marine Dataset, Strain Madness Dataset)
# Download GTDB genomics representative database for YACHT
cd $yacht_repo_loc
if [ ! -d $yacht_repo_loc/gtdb_reference ]; then
    mkdir $yacht_repo_loc/gtdb_reference;
    wget -P $yacht_repo_loc/gtdb_reference https://data.gtdb.ecogenomic.org/releases/release214/214.1/ar53_metadata_r214.tar.gz;
    wget -P $yacht_repo_loc/gtdb_reference https://data.gtdb.ecogenomic.org/releases/release214/214.1/bac120_metadata_r214.tar.gz;
    wget -P $yacht_repo_loc/gtdb_reference https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.zip;
fi

# Crate a reference dictionary matrix for YACHT
if [ ! -f $yacht_repo_loc/gtdb_ani_thresh_0.95_ref_matrix_processed.npz ]; then
    python $yacht_repo_loc/make_training_data_from_sketches.py --ref_file $yacht_repo_loc/gtdb_reference/gtdb-rs214-reps.k31.zip --ksize 31 --out_prefix 'gtdb_ani_thresh_0.95' --ani_thresh 0.95
fi

## run YACHT on CAMI2 data
# create a "results" directory for the output
if [ ! -d $benchmark_dir/results ]; then
    mkdir $benchmark_dir/results
fi
if [ ! -d $benchmark_dir/results/YACHT_results ]; then
    mkdir $benchmark_dir/results/YACHT_results
fi

# run YACHAT on Rhizosphere challenge data
if [ ! -d $benchmark_dir/results/YACHT_results/rhizosphere_data ]; then
    mkdir $benchmark_dir/results/YACHT_results/rhizosphere_data;
    ## create sketches of samples
    if [ ! -d $benchmark_dir/CAMI_data/rhizosphere_data/sketches ]; then
        mkdir $benchmark_dir/CAMI_data/rhizosphere_data/sketches;
        cd $benchmark_dir/CAMI_data/rhizosphere_data/sketches;
        parallel -j $cpu_num sourmash sketch dna -f -p k=31,scaled=1000,abund -o rhimgCAMI2_sample_{}.sig.zip $benchmark_dir/CAMI_data/rhizosphere_data/rhimgCAMI2_sample_{}.fq.gz ::: {0..20};
    fi
    ## run YACHT on the samples
    parallel -j $cpu_num python $yacht_repo_loc/run_YACHT.py --ref_matrix $yacht_repo_loc/gtdb_ani_thresh_0.95_ref_matrix_processed.npz --sample_file $benchmark_dir/CAMI_data/rhizosphere_data/sketches/rhimgCAMI2_sample_{}.sig.zip --ksize 31 --ani_thresh 0.95 --significance 0.99 --min_coverage 1 --outfile $benchmark_dir/results/YACHT_results/rhizosphere_data/rhimgCAMI2_sample_{}.csv ::: {0..20};
    ## convert the YACHT results to CAMI format
    parallel -j $cpu_num python $benchmark_dir/scripts/python_scripts/convert_to_CAMI_format.py --yacht_res $benchmark_dir/results/YACHT_results/rhizosphere_data/rhimgCAMI2_sample_{}.csv --gtdb_metadata_dir $yacht_repo_loc/gtdb_reference --outfile $benchmark_dir/results/YACHT_results/rhizosphere_data/rhimgCAMI2_sample_{}_cami_format.profile ::: {0..20};
fi

# run YACHT on Clinical pathogen detection challenge data
if [ ! -d $benchmark_dir/results/YACHT_results/pathogen_detection_data ]; then
    mkdir $benchmark_dir/results/YACHT_results/pathogen_detection_data;
    ## create sketches of samples
    if [ ! -d $benchmark_dir/CAMI_data/pathogen_detection_data/sketches ]; then
        mkdir $benchmark_dir/CAMI_data/pathogen_detection_data/sketches;
        cd $benchmark_dir/CAMI_data/pathogen_detection_data/sketches;
        sourmash sketch dna -f -p k=31,scaled=1000,abund -o patmgCAMI2.sig.zip $benchmark_dir/CAMI_data/pathogen_detection_data/*.fastq.gz;
    fi
    ## run YACHT on the samples
    python $yacht_repo_loc/run_YACHT.py --ref_matrix $yacht_repo_loc/gtdb_ani_thresh_0.95_ref_matrix_processed.npz --sample_file $benchmark_dir/CAMI_data/pathogen_detection_data/sketches/patmgCAMI2.sig.zip --ksize 31 --ani_thresh 0.95 --significance 0.99 --min_coverage 1 --outfile $benchmark_dir/results/YACHT_results/pathogen_detection_data/patmgCAMI2.csv;
    ## convert the YACHT results to CAMI format
    python $benchmark_dir/scripts/python_scripts/convert_to_CAMI_format.py --yacht_res $benchmark_dir/results/YACHT_results/pathogen_detection_data/patmgCAMI2.csv --gtdb_metadata_dir $yacht_repo_loc/gtdb_reference --outfile $benchmark_dir/results/YACHT_results/pathogen_detection_data/patmgCAMI2_cami_format.profile;
fi

# run YACHT on Challenge Marine Dataset
if [ ! -d $benchmark_dir/results/YACHT_results/marine_data ]; then
    mkdir $benchmark_dir/results/YACHT_results/marine_data;
    ## create sketches of samples
    if [ ! -d $benchmark_dir/CAMI_data/marine_data/sketches ]; then
        mkdir $benchmark_dir/CAMI_data/marine_data/sketches;
        cd $benchmark_dir/CAMI_data/marine_data/sketches;
        parallel -j $cpu_num sourmash sketch dna -f -p k=31,scaled=1000,abund -o marmgCAMI2_sample_{}.sig.zip $benchmark_dir/CAMI_data/marine_data/marmgCAMI2_sample_{}.fq.gz ::: {0..9};
    fi
    ## run YACHT on the samples
    parallel -j $cpu_num python $yacht_repo_loc/run_YACHT.py --ref_matrix $yacht_repo_loc/gtdb_ani_thresh_0.95_ref_matrix_processed.npz --sample_file $benchmark_dir/CAMI_data/marine_data/sketches/marmgCAMI2_sample_{}.sig.zip --ksize 31 --ani_thresh 0.95 --significance 0.99 --min_coverage 1 --outfile $benchmark_dir/results/YACHT_results/marine_data/marmgCAMI2_sample_{}.csv ::: {0..9};
    ## convert the YACHT results to CAMI format
    parallel -j $cpu_num python $benchmark_dir/scripts/python_scripts/convert_to_CAMI_format.py --yacht_res $benchmark_dir/results/YACHT_results/marine_data/marmgCAMI2_sample_{}.csv --gtdb_metadata_dir $yacht_repo_loc/gtdb_reference --outfile $benchmark_dir/results/YACHT_results/marine_data/marmgCAMI2_sample_{}_cami_format.profile ::: {0..9};
fi

# run YACHT on Strain Madness Dataset
if [ ! -d $benchmark_dir/results/YACHT_results/strain_madness_data ]; then
    mkdir $benchmark_dir/results/YACHT_results/strain_madness_data;
    ## create sketches of samples
    if [ ! -d $benchmark_dir/CAMI_data/strain_madness_data/sketches ]; then
        mkdir $benchmark_dir/CAMI_data/strain_madness_data/sketches;
        cd $benchmark_dir/CAMI_data/strain_madness_data/sketches;
        parallel -j $cpu_num sourmash sketch dna -f -p k=31,scaled=1000,abund -o strmgCAMI2_sample_{}.sig.zip $benchmark_dir/CAMI_data/strain_madness_data/strmgCAMI2_sample_{}.fq.gz ::: {0..99};
    fi
    ## run YACHT on the samples
    parallel -j $cpu_num python $yacht_repo_loc/run_YACHT.py --ref_matrix $yacht_repo_loc/gtdb_ani_thresh_0.95_ref_matrix_processed.npz --sample_file $benchmark_dir/CAMI_data/strain_madness_data/sketches/strmgCAMI2_sample_{}.sig.zip --ksize 31 --ani_thresh 0.95 --significance 0.99 --min_coverage 1 --outfile $benchmark_dir/results/YACHT_results/strain_madness_data/strmgCAMI2_sample_{}.csv ::: {0..99};
    ## convert the YACHT results to CAMI format
    parallel -j $cpu_num python $benchmark_dir/scripts/python_scripts/convert_to_CAMI_format.py --yacht_res $benchmark_dir/results/YACHT_results/strain_madness_data/strmgCAMI2_sample_{}.csv --gtdb_metadata_dir $yacht_repo_loc/gtdb_reference --outfile $benchmark_dir/results/YACHT_results/strain_madness_data/strmgCAMI2_sample_{}_cami_format.profile ::: {0..99};
fi
