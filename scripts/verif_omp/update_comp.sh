#!/bin/bash

# Compiler
COMP=$1

# Path to file.cfg
config_PATH="../../config/"
END="*.cfg"
CONFIG="/config/"

#LIST=${config_PATH}${COMP}${END} # without ${QP_ROOT}
LIST=${QP_ROOT}${CONFIG}${COMP}${END}

if [ -z "$1" ]
then
    echo "Give the compiler in argument"
else

    # List of the config files for the compiler
    #list_files=$(ls ../../config/$comp*.cfg) #does not give the right list
    list_files=${LIST}
    echo "Files that will be modified:"
    echo $list_files
   
    # Flags that must be added 
    FLAGS=$(./check_required_setup.sh $COMP)

    # Add the flags
    for file in $list_files
    do
        echo $file
        BASE="IRPF90_FLAGS : --ninja"
        ACTUAL=$(grep "$BASE" $file)
        
        # To have only one time each flag	
	grep " -DSET_MAX_ACT" $file && ${ACTUAL/" -DSET_MAX"/""}
	grep " -DSET_NESTED" $file && ${ACTUAL/" -DSET_NESTED"/""}
        SPACE=" "

	NEW=${ACTUAL}${SPACE}${FLAGS}

	# Debug
        #echo ${NEW}

        sed "s/${ACTUAL}/${NEW}/" $file
        # -i # to change the files
    done

fi
