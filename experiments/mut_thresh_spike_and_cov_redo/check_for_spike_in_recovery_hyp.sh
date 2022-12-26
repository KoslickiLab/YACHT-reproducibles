#!/usr/bin/env bash
#set -e
#set -u
set -o pipefail

# For each of the experiments, find who they are supposed to be similar to and check if it's in the sample
# prep the output csv
echo "spike_coverage@coverage_threshold@max_ani@rel_ab@mut_thresh@spike_md5short@spike_name@match_name" > results.csv
coverageValues=("1" "0.1" "0.01" "0.001")
#mutThreshs=("0.1" "0.05" "0.01" "0.001")
mutThreshes=("0.1")
for mutThresh in "${mutThreshes[@]}"
do
	for spikeCov in "${coverageValues[@]}"
	do
        	for covThresh in "${coverageValues[@]}"
        	do
        	        numLines=$(wc -l sigs_md5_to_accession_to_gtdb_location_mut_${mutThresh}.txt | cut -d' ' -f1)
        	        it=0
                	for line in `tail -n +2 sigs_md5_to_accession_to_gtdb_location_mut_${mutThresh}.txt`
                	do
                	        #increment the iterator
                	        it=$((it+1))
                	        # if the iterator is a multiple of 100, print a message
                	        if [ $((it%100)) -eq 0 ]; then
                	                echo "On line ${it} of ${numLines}"
                	        fi
                        	md5short=$(echo ${line} | cut -d',' -f2)
                          md5long=$(echo ${line} | cut -d',' -f1)
                          match=$(grep -P -m 1 ${md5long} in_gtdb_similar_to_EU_not_in_sample_mut_${mutThresh}_@.tsv)
                          # Just take the accession, so we don't need to deal with commas in names
			                    gtdbName=$(echo $match | cut -d'@' -f1 | cut -d' ' -f1)
                          gtdbMD5=$(echo $match | cut -d'@' -f2)
                          maxANI=$(echo $match | cut -d'@' -f3)
                          matchingOrg=$(echo $match | cut -d'@' -f4)
                          #read -d'@' -r gtdbName gtdbMD5 maxANI matchingOrg <<< $(grep -m 1 ${md5long} in_gtdb_similar_to_EU_not_in_sample_@.tsv)
                          #echo $gtdbName
                          #echo $gtdbMD5
                          #echo $maxANI
                          #echo $matchingOrg
                          #exit
                          gtdbShortName=$(echo ${gtdbName} | cut -d'.' -f1)
                          #head -n 1 EU_on_spikes_cov_${spikeCov}/${md5short}_${spikeCov}X_cov_thresh_${covThresh}X.csv
                          #grep ${matchingOrg} EU_on_spikes_cov_${spikeCov}/${md5short}_${spikeCov}X_cov_thresh_${covThresh}X.csv
                          relAb=$(grep ${matchingOrg} EU_on_spikes_cov_${spikeCov}_hyp/${md5short}_${spikeCov}X_cov_thresh_${covThresh}X_mut_thresh_${mutThresh}_hyp.csv | cut -d',' -f 13)
                          if [ $? -eq 0 ]; then
                            #echo ${relAb}
                            true
                          else
                                relAb="NaN"
                          fi
                          echo "${spikeCov}@${covThresh}@${maxANI}@${relAb}@${mutThresh}@${md5short}@${gtdbName}@${matchingOrg}" >> results.csv
			            done
		      done
	 done
done