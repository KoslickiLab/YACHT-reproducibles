# YACHT (Benchmarking and Proof-of-Concept Experiments)

**_PLEASE NOTE: this repo is created for the proof-of-concept and benchmarking experiments of YACHT. For updates on the production-level code, please follow [this repo](https://github.com/KoslickiLab/YACHT/)_**

## Installation
Please install the necessary packages via [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) following the commands below:

```bash
# Clone the repo
git clone https://github.com/KoslickiLab/YACHT.git
cd YACHT

# Create a conda environment for YACHT
conda env create -f env/yacht_env.yaml
conda activate yacht
```

## Benchmarking Experiments

To evaluate the performance of YACHT, we utilize the robust [public datasets](https://data.cami-challenge.org/participate) from the [Critical Assessment of Metagenome Interpretation II (CAMI II)](https://www.nature.com/articles/s41592-022-01431-4) for the tasks of taxonomic profiling and clinical pathogen detection, which contains the rhizosphere data, marine data, strain madness data, and clinical pathogen data. These datasets include the rhizosphere data, marine data, strain madness data, and clinical pathogen data, providing high-quality benchmarks for evaluating metagenomic analysis tools like YACHT. We also leverage the CAMI-official profiling assessment tool [OPAL](https://github.com/CAMI-challenge/OPAL/tree/master) to compare YACHT's performance against the state-of-the-art (SOTA) tools (e.g., Bracken, Metaglin, mOTUs, MetaPhlAn, CCMetagen, NBC++, MetaPhyler, LSHVec) suggested by CAMI II. 

### How to reproduce the evaluation results
We provide two bash scripts under `./benchmark/scripts` folder to reproduce our evaluation results. To run these scripts, please make sure your have installed the Conda environnment suggested above. After that, run the following instruction:
1. Git clone the [production-level YACHT repositoary](https://github.com/KoslickiLab/YACHT):
```bash
git clone https://github.com/KoslickiLab/YACHT.git
```

2. Download CAMI2 datasets 
```bash
# download_cami2_data.sh <benchmark_dir> <cpu_num>
bash ./benchmark/scripts/bash_scripts/download_cami2_data.sh <path_to_YACHT_experiment>/benchmark 50
```

3. Run YACHT on the CAMI datasets
```bash
# run_YACHT.sh <yacht_repo_loc> <benchmark_dir> <cpu_num>
bash ./benchmark/scripts/bash_scripts/run_YACHT.sh <path_to_YACHT_experiment>/YACHT <path_to_YACHT_experiment>/benchmark 20
```

4. Run OPAL on the YACHT results
```bash
# git clone OPAL 
git clone https://github.com/CAMI-challenge/OPAL
# run_OPAL.sh <opal_repo_loc> <benchmark_dir> <cpu_num>
bash ./benchmark/scripts/bash_scripts/run_OPAL.sh <path_to_YACHT_experiment>/OPAL <path_to_YACHT_experiment>/benchmark 20
```

## Proof-of-Concept Experiments

### Creating a reference dictionary matrix (`ref_matrix.py`):
```bash 
python ref_matrix.py --ref_file '../ForSteve/ref_gtdb-rs207.genomic-reps.dna.k31.zip' --out_prefix 'test2_' --N 20
```

### Computing relative abundance of organisms (`recover_abundance.py`):
```bash
python recover_abundance.py --ref_file 'test2_ref_matrix_processed.npz' --sample_file '../ForSteve/sample.sig' --hash_file 'test2_hash_to_col_idx.csv' --org_file 'test2_processed_org_idx.csv' --w 0.01 --outfile 'test2_recovered_abundance.csv'
```

### Basic workflow
1. run ```python ref_matrix.py --ref_file 'tests/testdata/20_genomes_sketches.zip' --out_prefix 'tests/unittest_'``` . This 
should generate 4 files in the tests folder.
2. run ```python recover_abundance.py --ref_file 'tests/unittest_ref_matrix_processed.npz' --sample_file 
   'tests/testdata/sample.sig' --hash_file 'tests/unittest_hash_to_col_idx.csv' --org_file 'tests/unittest_processed_org_idx.csv' --w 0.01 --outfile 'tests/unittest_recovered_abundance.csv'``` . Should create a file `tests/unittest_recovered_abundance.csv` which should be all zeros.
3. run the same command as above, but with `--w 0.0001`. Should overwrite `tests/unittest_recovered_abundance.csv` with a 
   6 in the 19th row
