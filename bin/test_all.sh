#!/bin/bash
OUT="/tmp/scripts.out"
#nosetests scripts > ${OUT} 2>&1
grep "Exception" ${OUT} | grep -v "CVODE"
grep "Error" ${OUT} | grep -v "CVODE"
grep "Fail" ${OUT}

