#!/bin/bash

# Project directory naming #

# `cascade_version` (which network and corresponding directory to use): {1.0, 2.0}
# `train` (steady state train vs proliferation train): {ss, rand}
# `sim_num` (number of simulations to run): {AnyListOfNumbers > 0}
# `attr_tool` (attractor-tool): {fixpoints = biolqm_stable_states, bnet = bnet_reduction_reduced, traps = biolqm_trapspaces}
# `synergy_method` (Drabme's configuration option): {hsa, bliss}

for cascade_version in 1.0 2.0
do
    for train in ss rand
    do
	if [ $train == "ss" ]
       	then
	    cat ags_cascade_$cascade_version/steadystate > ags_cascade_$cascade_version/training
	elif [ $train == "rand" ]
	then
	    cat ags_cascade_$cascade_version/random_train > ags_cascade_$cascade_version/training
        fi
	for sim_num in 50 150
	do
	    sed_str='s/simulations:.*/simulations:\t'${sim_num}'/'
	    sed -i $sed_str ags_cascade_$cascade_version/config
	    for attr_tool in fixpoints
	    do
	        if [ $attr_tool == "bnet" ]
		then
		    sed -i 's/attractor_tool:.*/attractor_tool:\tbnet_reduction_reduced/' ags_cascade_$cascade_version/config
		elif [ $attr_tool == "fixpoints" ]
                then
                    sed -i 's/attractor_tool:.*/attractor_tool:\tbiolqm_stable_states/' ags_cascade_$cascade_version/config
	        elif [ $attr_tool == "traps" ]
		then	
		    sed -i 's/attractor_tool:.*/attractor_tool:\tbiolqm_trapspaces/' ags_cascade_$cascade_version/config
                fi

		for synergy_method in hsa bliss
		do
		  sed_str='s/synergy_method:.*/synergy_method:\t'${synergy_method}'/'
		  sed -i $sed_str ags_cascade_$cascade_version/config
		
	          # Run Launcher
		  project_name=cascade_${cascade_version}_${train}_${sim_num}sim_${attr_tool}_${synergy_method}
		  echo START: $project_name
		  start=`date +%s`
		  java -cp target/synergy-1.2.1-jar-with-dependencies.jar eu.druglogics.synergy.Launcher --inputDir=ags_cascade_${cascade_version} --project=$project_name > /dev/null 2>&1
		  runtime=$(($(date +%s)-$start))
		  echo -e Execution Time: "$(($runtime / 60)) minutes\n"
		done
	    done
        done
    done
done

