#!/bin/bash
# Creates a script from a ipython notebook
# input: file.ipynb
# output: file.py
INPUT=$1
echo "*** mk_script for ${INPUT}"
FILE=${INPUT/.ipynb/}
echo $FILE
NOTEBOOK="$FILE.ipynb"
PYTHON="${FILE}.py"
TEXT="${FILE}.txt"
TMP_PYTHON="/tmp/${FILE}.py"
DEP="/tmp/${FILE}"
jupyter nbconvert --to script "${NOTEBOOK}"
if [ ! -f ${PYTHON} ]
    then
        mv ${TEXT} ${PYTHON}
fi
echo "No dependences" > ${DEP}
grep "pip install" ${PYTHON} > ${DEP}
grep -v "pip install" ${PYTHON} > ${TMP_PYTHON}
cp ${TMP_PYTHON} ${PYTHON}
grep -v "get_ipython" ${PYTHON} > ${TMP_PYTHON}
cp ${TMP_PYTHON} ${PYTHON}
sed 's/IS_COLAB = True/IS_COLAB = False/' ${PYTHON} > ${TMP_PYTHON}
cp ${TMP_PYTHON} ${PYTHON}
echo ""
echo "**Below are the dependencies for this ${PYTHON}"
cat ${DEP}
echo ""
