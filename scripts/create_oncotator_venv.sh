#!/bin/bash
# Adapted from internal Broad cgap project.
# Have script stop if there is an error
set -e

#################################################
# Parsing arguments
#################################################

while getopts ":e:ckstm" option; do
	case "$option" in
		e) ENV="$OPTARG" ;;
		c) FLAGS="archflags" ;;
		k) COMPIL="skip" ;;
		s) PYVCF="skip" ;;
		t) TRAVIS=true ;;
		m) MAC=true ;;
	esac
done

if [ -z "$ENV" ] && [ ! $TRAVIS ]; then
	printf "Option -e requires an argument.\n \
Usage: %s: -e <ENV> [-cks] \n \
where <ENV> is the name to associate with this environment \n \
(e.g. create_oncotator_venv.sh oncotator_test_env)\n \
Optional arguments:  \n \
-c \t trigger the workaround for the XCode 5.1.1 compilation bug \n \
   \t (see documentation for details)  \n \
-k \t skip packages that require compilation  \n \
-s \t skip packages that are not available through PyPI \n \
-t \t run in travis installation mode \n \
-m \t create venv on Mac (i.e. skip ngslib installation) \n" $0

	exit 1
fi

if [ ! -z "$FLAGS" ]; then
	printf "Option -c specified -- the ARCHFLAGS workaround will be applied.\n"
fi

if [ ! -z "$COMPIL" ]; then
	printf "Option -k specified -- packages that require compilation will be skipped.\n"
fi

if [ ! -z "$PYVCF" ]; then
	printf "Option -s specified -- packages that cannot be obtained through PyPI will be skipped.\n"
fi

if [ ! -z "$TRAVIS" ]; then
	printf "TRAVIS environment variable is set.  Do not activate/deactivate virtual envs.\n"
fi

SKIP_MSG="Skipping... Make sure to install these packages manually after the script has finished. "

#################################################
# Create the V-ENV
#################################################

if [ ! $TRAVIS ]; then
	# Create and activate a test environment
	virtualenv $ENV
	source $ENV/bin/activate

	echo " "
	echo "Virtual environment created and activated in $ENV."
	echo "Now attempting to install packages into the virtual environment."
else
	which python
	python --version
fi

echo "Update Pip"
pip --version
pip install -U pip
pip --version

#################################################
# Installations that require compilation
#################################################

if [ "$COMPIL" == "skip" ];
then
	echo $SKIP_MSG
else
	echo "Attempting to install packages that require compilation. If this fails, try again with the flag -c added to the script command. If that still does not work, you will need to install them manually."

	for C_PACKAGE in biopython cython numpy pandas sqlalchemy
	do
		echo " "
		echo "$C_PACKAGE =========================="
		if [ "$FLAGS" == "archflags" ]; then
			env ARCHFLAGS="-Wno-error=unused-command-line-argument-hard-error-in-future" pip install $C_PACKAGE
		else
			pip install --no-use-wheel $C_PACKAGE
		fi
		echo "OK"
	done

	if [ ! $MAC ]; then
		echo " "
		echo "ngslib =========================="
		if [ "$FLAGS" == "archflags" ]; then
			env ARCHFLAGS="-Wno-error=unused-command-line-argument-hard-error-in-future" pip install ngslib
		else
			pip install --no-use-wheel ngslib
		fi
		echo "OK"
	fi		
fi

#################################################
# Easy installations
#################################################

echo " "
echo "Installing dependencies that can be obtained from pypi"

for PACKAGE in bcbio-gff nose shove python-memcached natsort more-itertools enum34
do
	echo " "
	echo "$PACKAGE =========================="
	pip install -U --no-use-wheel $PACKAGE
	echo "OK"
done


#################################################
# Tricky installations
#################################################

echo " "
echo "pyvcf ========================="

if [ "$PYVCF" == "skip" ];
then
	echo $SKIP_MSG
else
	echo "Attempting to install a package that cannot be obtained from PyPI. If this fails, you will need to install it manually after the script has run. "
	# This one is a little complicated
	echo "Retrieving mgupta (aka elephanthunter) fork of PyVCF"
	wget --no-check-certificate 'https://github.com/elephanthunter/PyVCF/archive/master.zip'

	if [ -f "master" ];
	then
	   mv master master.zip
	else
	   echo "No master found, assuming master.zip"
	fi

	unzip master.zip && cd PyVCF-master && python setup.py install && cd .. && rm -Rf PyVCF-master && rm -f master.* && rm -f master*

	echo "OK."
fi

#################################################
# All done!
#################################################

echo "NOTE: Oncotator has not been installed, only the dependencies. You MUST still install Oncotator manually. "
if [ ! $TRAVIS ]; then
	echo "Deactivating"
	deactivate
fi
