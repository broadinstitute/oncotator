#!/bin/bash
# Adapted from internal Broad cgap project.
# Have script stop if there is an error
set -e

if [ -z "$1" ]
  then
    echo "usage: create_cgap_venv.sh <ENV>";
    echo "    where <ENV> is the name to associate with this environment (e.g. create_oncotator_venv.sh oncotator_test_env)";
    exit 1
fi



# Create and activate a test environment
virtualenv $1
source $1/bin/activate

#################################################
# Manual installations (third party automatic download in pypi is broken)
#################################################
echo "Attempting to install packages that cannot be installed properly from pypi.  Any other oncotator dependencies can be installed through installing oncotator or manually"
echo " "
echo "nose =========================="
pip install nose

echo " "
echo "numpy ========================="
pip install numpy

echo " "
echo "PyVCF ========================="
# This one is a little more complicated
echo "Retrieving mgupta (aka elephanthunter) fork of PyVCF"
wget --no-check-certificate 'https://github.com/elephanthunter/PyVCF/archive/master.zip'
#mv master master.zip
unzip master.zip && cd PyVCF-master && python setup.py install && cd .. && rm -Rf PyVCF-master && rm -f master.* && rm -f master*

echo " "
echo "BioPython  ========================="
pip install biopython

echo " "
echo "Installing most dependencies that can be downloaded from pypi"

echo "bcbio-gff ========================="
pip install bcbio-gff
echo " "
echo "pysam ========================="
pip install pysam
echo " "
echo "pandas ========================="
pip install pandas
echo " "
echo "cython ========================="
pip install cython
echo " "
echo "shove ========================="
pip install shove
echo " "
echo "sqlalchemy ========================="
pip install sqlalchemy
echo " "
echo "python-memcached ========================="
pip install python-memcached

echo " NOTE:  Some dependencies may not have been installed, you should still install oncotator normally."

deactivate
