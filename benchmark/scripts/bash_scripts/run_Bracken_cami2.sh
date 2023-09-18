#! /bin/bash
# run pathogen detection experiments based on CAMI2 data with Bracken tool

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: run_Bracken_cami2.sh <yacht_reproducibles_dir> <cpu_num>"
    exit 1
fi
yacht_reproducibles_dir=$1
benchmark_dir=${yacht_reproducibles_dir}/benchmark
cpu_num=$2
dbname=standard_db
KRAKEN_DB=${yacht_reproducibles_dir}/kraken2/${dbname}

# ## Download and install Bracken
# echo "Downloading and installing Bracken..."
# cd $yacht_reproducibles_dir
# git clone https://github.com/jenniferlu717/Bracken.git
# cd Bracken
# bash install_bracken.sh

# ## Download and install Kraken2
# echo "Downloading and installing Kraken2..."
# cd $yacht_reproducibles_dir
# git clone https://github.com/DerrickWood/kraken2.git
# cd kraken2
# ./install_kraken2.sh .

## Creat a result folder
echo "Creating a result folder..."
if [ ! -d "${benchmark_dir}/results/bracken_result" ]; then
    mkdir ${benchmark_dir}/results/bracken_result
fi

# ## Step 0: Build a Kraken2 database
# echo "Building a Kraken2 database..."
# ${yacht_reproducibles_dir}/kraken2/kraken2-build --standard --threads $cpu_num --db $dbname
# # ./${working_dir}/kraken2/kraken2-build --download-taxonomy --db 'customized_db'
# # ./${working_dir}/kraken2/kraken2-build --download-library archaea --db "customized_db"
# # ./${working_dir}/kraken2/kraken2-build --download-library fungi --db "customized_db"
# # ./${working_dir}/kraken2/kraken2-build --download-library bacteria --db "customized_db"
# # ./${working_dir}/kraken2/kraken2-build --download-library viral --db "customized_db"

# ## Step 1: Generate the Bracken database file
# echo "Generating the Bracken database file..."
# ${yacht_reproducibles_dir}/Bracken/bracken-build -d ${KRAKEN_DB} -t $cpu_num -x ${yacht_reproducibles_dir}/kraken2

# ## Step 2: Run Kraken 2 AND Generate a report file
# echo "Running Kraken 2 AND Generating a report file..."
# ${yacht_reproducibles_dir}/kraken2/kraken2 --db ${KRAKEN_DB} --threads $cpu_num --report ${benchmark_dir}/results/bracken_result/k2_report.txt --output ${benchmark_dir}/results/bracken_result/k2_output.txt --gzip-compressed --paired ${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R1.fastq.gz ${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R2.fastq.gz 
# ${yacht_reproducibles_dir}/Bracken/bracken -d ${KRAKEN_DB} -i ${benchmark_dir}/results/bracken_result/k2_report.txt -o ${benchmark_dir}/results/bracken_result/bracken_output.txt -t $cpu_num -w ${benchmark_dir}/results/bracken_result/bracken_report.txt

${yacht_reproducibles_dir}/kraken2/kraken2 --db ${KRAKEN_DB} --threads $cpu_num --report ${benchmark_dir}/results/bracken_result/k2_report.txt --output ${benchmark_dir}/results/bracken_result/k2_output.txt --gzip-compressed /home/grads/cqm5886/interm_files/YACHT-reproducibles/real_world_experiment/data/reads_fastq/ERR4402517/ERR4402517.fastq.gz 
${yacht_reproducibles_dir}/Bracken/bracken -d ${KRAKEN_DB} -i ${benchmark_dir}/results/bracken_result/k2_report.txt -o ${benchmark_dir}/results/bracken_result/bracken_output.txt -t $cpu_num -w ${benchmark_dir}/results/bracken_result/bracken_report.txt

