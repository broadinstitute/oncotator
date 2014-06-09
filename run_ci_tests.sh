#!/bin/bash
# Have script stop if there is an error
echo "Starting CI testing..."
set -e

# Make sure that a db-dir is specified, if not, set it to the default for the oncotator CI server at Broad.
if [ -z "$1" ]
  then
    echo "usage: run_ci_tests.sh <DB_DIR>";
    echo "    where <DB_DIR> is the location of your datasources ";
    exit 1
  else
    DB_DIR=${1}
fi


echo "This script must be run from the same directory as setup.py"

#echo "PATH=/xchip/tcga/Tools/oncotator/python_2.7.6_May192014/bin/:$PATH" > tmp_rc
#source tmp_rc
which python
python --version

VENV=oncotator_nosetest_env_auto
rm -Rf $VENV
bash ./scripts/create_oncotator_venv.sh -e $VENV
source $VENV/bin/activate
echo "VENV python:"
which python


# nosetests will return a non-zero error code if any unit test fails, so we do not want to stop this script
#   So we remove the error catching temporarily
python setup.py install

mkdir -p out

# Create a proper config file
sed -r "s:dbDir=MY_DB_DIR:dbDir=${DB_DIR}:g" test/configs/personal-test.config.template >test/configs/personal-test.config

# Attempt simple command line functionality
echo "Attempting a few command line calls (oncotator only) to make sure there are no egregious errors in the Oncotator CLI."
echo "== Just running --help ==="
oncotator --help

echo "== Indel maflite 2 tcga maf test ==="
oncotator -v --no-multicore --db-dir=${DB_DIR} test/testdata/maflite/Patient0.indel.maf.txt $VENV/test_Patient0.indel.maf.txt hg19

echo "== SNV maflite 2 tcga maf test ==="
oncotator -v --no-multicore --db-dir=${DB_DIR} test/testdata/maflite/Patient0.snp.maf.txt $VENV/test_Patient0.snp.maf.txt hg19

echo "== maflite 2 vcf infer genotypes test"
oncotator -v --no-multicore --db-dir=${DB_DIR} --infer_genotypes=true --input_format=MAFLITE --output_format=VCF test/testdata/maflite/Patient0.snp.maf.txt $VENV/test_Patient0.snp.maf.vcf hg19

########

set +e
# Do not use multiprocess mode with profiling or coverage.  Bug in nosetests also disallows --processes and --with-xunit
#  --processes=4 --process-timeout=480  --process-restartworker
nosetests --all-modules --exe --with-xunit --xunit-file=${PWD}/nosetests.xml -w test -v
set -e

echo "Deactivating and deleting test python virtual environment"
deactivate
rm -Rf $VENV
