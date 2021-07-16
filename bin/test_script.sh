#!/bin/bash
# Tests scripts in a directory
DIR=`pwd`
echo ""
echo "***Processing ${DIR}"
for f in `ls *.py`
    do
        OUT="/tmp/${f/.py/.out}"
        echo "   Processing $f into ${OUT}."
        python $f > ${OUT} 2>&1
    done
