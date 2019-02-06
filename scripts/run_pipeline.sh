#!/bin/bash

while read parameter ; do

	if [[ ! "$parameter" =~ "platform" ]]; then

		echo Selected parameters are : $parameter;
		platforms=$(echo $parameter | awk '{print $1}')
		lengths=$(echo $parameter | awk '{print $2}')
		coverages=$(echo $parameter | awk '{print $3}')
		meansizes=$(echo $parameter | awk '{print $4}')
		stdevs=$(echo $parameter | awk '{print $5}')
		minqualities=$(echo $parameter | awk '{print $6}')
		maxinsertions=$(echo $parameter | awk '{print $7}')


		# now run for each of these parameters

		for platform in $(echo $platforms | sed 's/,/ /g'); do
			for length in $(echo $lengths | sed 's/,/ /g'); do
				for coverage in $(echo $coverages | sed 's/,/ /g'); do
					for meansize in $(echo $meansizes | sed 's/,/ /g'); do
						for stdev in $(echo $stdevs | sed 's/,/ /g'); do
							for minquality in $(echo $minqualities | sed 's/,/ /g'); do
								for maxinsertion in $(echo $maxinsertions | sed 's/,/ /g'); do

									datafolder=/tsl/scratch/shrestha/assembly_comparison/results/illumina/platform${platform}/length${length}/coverage${coverage}/meansize${meansize}/stdev${stdev}/minquality${minquality}/maxinsertion${maxinsertion}
									mkdir -p $datafolder

									#simulate data
									echo "Simulating illumina data with these parameters platform${platform}_length${length}_coverage${coverage}_meansize${meansize}_stdev${stdev}_minquality${minquality}_maxinsertion${maxinsertion} "

									# prefix for output simulation data
									prefix=/tsl/scratch/shrestha/assembly_comparison/results/illumina_simulation_data/simulation_of_platform${platform}_length${length}_coverage${coverage}_meansize${meansize}_stdev${stdev}_minquality${minquality}_maxinsertion${maxinsertion}_

									if test -e ${prefix}complete && test ${prefix}complete -nt ${prefix}1.fq; then
										echo "simulation of ${prefix}1.fq and it's paired end reads completed";
									else

										echo Simulating data with prefix ${prefix}
										bash ~/workarea/assembly_comparison/scripts/simulate_art_illumina_reads.sh /hpc-home/shrestha/workarea/assembly_comparison/data/PGSC_DM_v3_scaffolds.fasta $platform $length $coverage $meansize $stdev $minquality $maxinsertion $prefix

										# create a file that marks file simulation was completed.
										touch ${prefix}complete
										
										# sleep for a while now 
										echo "Sleeping 2hrs. current time is " $(date)
										sleep 3600

									fi

									# do assembly and other analyses
									filecounter=0

									cd $datafolder
									masterscript=${datafolder}/analysis_masterscript.sh

									#if test -e $masterscript; then echo "$masterscript exists"; continue; fi

									for kmer in `seq 21 4 55`; do
										for mincoverage in 3 5 10; do
											for minkmercoverage in 10 20 30 40 50; do
												filecounter=$((filecounter + 1))
												assemblyfolder=${datafolder}/kmer${kmer}/mincoverage${mincoverage}/minkmercoverage${minkmercoverage}
												mkdir -p $assemblyfolder
												cd $assemblyfolder

												echo '#!/bin/bash' > ${datafolder}/analysis_${filecounter}.sh
												echo -e "cd ${assemblyfolder};\n\nif test -e analysis_complete; then\necho Analysis already completed; exit;\nelse\nbash ~/workarea/assembly_comparison/scripts/run_abyss_blast_quast_transrate.sh $kmer $mincoverage $minkmercoverage ${prefix}1.fq ${prefix}2.fq > ${assemblyfolder}/analysis.log && touch analysis_complete;\nfi "  >> ${datafolder}/analysis_${filecounter}.sh
											done
										done
									done
									cd $datafolder
									masterscript=${datafolder}/analysis_masterscript.sh

									# submit job array
									echo '#!/bin/bash' > ${masterscript}
									echo "srun bash ${datafolder}/analysis_\${SLURM_ARRAY_TASK_ID}.sh > ${datafolder}/analysis_output_\${SLURM_ARRAY_JOB_ID}.log" >> ${masterscript}
									echo "sbatch the job array now. Master script is: $masterscript"
									sbatch  -o ${assemblyfolder}/masterscript_analysis.log  --array=1-${filecounter} --mem 60G --cpus 2 -J kmer${kmer}_mincoverage${mincoverage}_minkmercoverage${minkmercoverage} ${masterscript}


								done
							done
						done
					done
				done
			done
		done

	fi;



done < ~/workarea/assembly_comparison/data/parameters_for_data_simulation.txt
