#! /bin/bash
# run pathogen detection experiments based on CAMI2 data with MetaPhAn tool

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: run_MetaPhAn_cami2.sh <yacht_reproducibles_dir> <cpu_num>"
    exit 1
fi
yacht_reproducibles_dir=$1
benchmark_dir=${yacht_reproducibles_dir}/benchmark
cpu_num=$2
## check for the latest version of MetaPhlAn database for mpa v3 in http://cmprod1.cibio.unitn.it/biobakery3/metaphlan_databases/ because we need to run with "--add_viruses --mpa3"
metaphlan_db_version="mpa_v31_CHOCOPhlAn_201901"

## Creat a result folder
echo "Creating a result folder..."
if [ ! -d "${benchmark_dir}/results/metaphlan_result" ]; then
    mkdir ${benchmark_dir}/results/metaphlan_result
fi

## Run MetaPhlAn
echo "Unzipping the input file..."
if [ ! -e "${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R1.fastq" ]; then
    gunzip -c ${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R1.fastq.gz > ${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R1.fastq
    metagenome_1=${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R1.fastq
else
    metagenome_1=${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R1.fastq
fi
if [ ! -e "${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R2.fastq" ]; then
    gunzip -c ${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R2.fastq.gz > ${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R2.fastq
    metagenome_2=${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R2.fastq
else
    metagenome_2=${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R2.fastq
fi
# combine two fastq files into one
echo "Combining two fastq files into one..."
if [ ! -e "${benchmark_dir}/CAMI_data/pathogen_detection_data/combined.fastq" ]; then
    cat ${metagenome_1} ${metagenome_2} > ${benchmark_dir}/CAMI_data/pathogen_detection_data/combined.fastq
    metaphlan_in=${benchmark_dir}/CAMI_data/pathogen_detection_data/combined.fastq
else
    metaphlan_in=${benchmark_dir}/CAMI_data/pathogen_detection_data/combined.fastq
fi


echo "Downloading MetaPhlAn database..."
if [ ! -d "${benchmark_dir}/results/metaphlan_result/metaphlan_databases" ]; then
    mkdir -p ${benchmark_dir}/results/metaphlan_result/metaphlan_databases
    wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/${metaphlan_db_version}.tar -P ${benchmark_dir}/results/metaphlan_result/metaphlan_databases
    cd ${benchmark_dir}/results/metaphlan_result/metaphlan_databases
    tar -xvf ${metaphlan_db_version}.tar
    ls *.bz2 | xargs -I {} bzip2 -d {}
    ls *.fna | parallel -j $cpu_num "outfile=\$(echo {} | sed 's/.fna/.bwt2/g'); bowtie2-build {} \$outfile --large-index --threads \$cpu_num"
    rm ${metaphlan_db_version}.tar
    cp ${metaphlan_db_version}.pkl ${metaphlan_db_version}.bwt2.pkl
fi

echo "Running MetaPhlAn..."
metaphlan ${metaphlan_in} --input_type fastq --add_viruses --mpa3 --nproc $cpu_num --CAMI_format_output --index ${metaphlan_db_version}.bwt2 --bowtie2db ${benchmark_dir}/results/metaphlan_result/metaphlan_databases -o ${benchmark_dir}/results/metaphlan_result/metaphlan_output.tsv
