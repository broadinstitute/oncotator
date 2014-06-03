#!/bin/bash
# Adapted from internal Broad cgap project.
# Have script stop if there is an error
set -e

#################################################
# Parsing arguments
#################################################

while getopts ":e:cks" option; do
	case "$option" in 
		e) ENV="$OPTARG" ;;
		c) FLAGS="archflags" ;;
		k) COMPIL="skip" ;;
		s) PYVCF="skip" ;;	
	esac 
done

if [ -z "$ENV" ]; then
	printf "Option -e requires an argument.\n \
Usage: %s: -e <ENV> [-cks] \n \
where <ENV> is the name to associate with this environment \n \
(e.g. create_oncotator_venv.sh oncotator_test_env)\n \
Optional arguments:  \n \
-c \t trigger the workaround for the XCode 5.1.1 compilation bug \n \
   \t (see documentation for details)  \n \
-k \t skip packages that require compilation  \n \
-s \t skip packages that are not available through PyPI \n" $0
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
 
SKIP_MSG="Skipping... Make sure to install these packages manually after the script has finished. "

#################################################
# Create the V-ENV
#################################################

# Create and activate a test environment
virtualenv $ENV
source $ENV/bin/activate

echo " "
echo "Virtual environment created and activated in $ENV."
echo "Now attempting to install packages into the virtual environment."
which python
python --version


#################################################
# Easy installations
#################################################

echo " "
echo "Installing dependencies that can be obtained from pypi"

for PACKAGE in bcbio-gff nose shove python-memcached natsort
do 
	echo " "
	echo "$PACKAGE =========================="
	pip install -U $PACKAGE
	echo "OK"
done

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
			pip install $C_PACKAGE
		fi
		echo "OK"
	done
	
	echo " "
	echo "pysam ========================="
	if [ "$FLAGS" == "archflags" ]; then
		env ARCHFLAGS="-Wno-error=unused-command-line-argument-hard-error-in-future" pip install -I --allow-unverified pysam pysam==0.7.5
	else 
		pip install -I pysam==0.7.5
	fi
	echo "OK"
fi

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

deactivate
