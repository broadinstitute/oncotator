#!/bin/bash

# This script will CHECK all files in the Oncotator for the license information
# it is meant to be run by bamboo as a sanity check that our repo contains the
# correct license information in all files.
#
# script must be run from the $GIT_ROOT
#
# Adapted from a script for the GATK:
# author: Mauricio Carneiro
# date: 1/9/13

echo "Checking all licenses in Oncotator... ";

result=0
ls oncotator/**/*.py     | python private/python/licensing/CheckLicense.py
if [[ "$?" == "255" ]]
then
    result=1
fi

ls scripts/**/*.py     | python private/python/licensing/CheckLicense.py
if [[ "$?" == "255" ]]
then
    result=1
fi

ls test/**/*.py  | python private/python/licensing/CheckLicense.py
if [[ "$?" == "255" ]]
then
    result=1
fi

exit $result

