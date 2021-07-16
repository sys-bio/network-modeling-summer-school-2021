#!/bin/bash
# Builds all of the scripts
# Names of files and directories should not contain spaces
read -p "*** This script clobbers the scripts/ direction. Hit enter to continue"

START_DIR=`pwd`
BASE_DIR="${START_DIR}/scripts"  # Target directory
CMD="${START_DIR}/mk_script.sh"

rm -rf scripts
cp -r notebooks scripts

# Creates a script file
function create() {
    if [ -z "$1" ]
        then
            DIR="${BASE_DIR}"
        else
            OLD_DIR="${BASE_DIR}/$1"
            DIR="${BASE_DIR}/test$1"
            mv "${OLD_DIR}" "${DIR}"
    fi
    echo "Processing ${DIR}"
    cd "${DIR}"
    for f in `ls *.ipynb`
      do
        bash "${START_DIR}/mk_script.sh" "$f"
        rm $f
        echo $f
      done
    for f in *.py
      do
        mv $f test$f
      done
}

# Top directory
create
cd ${BASE_DIR}
LIST=`echo */ | grep -v ".py"`
for d in ${LIST}
    do
        create $d
    done
