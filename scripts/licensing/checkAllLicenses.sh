#!/bin/zsh

# This script will CHECK all files in the GATK for the license information
# it is meant to be run by bamboo as a sanity check that our repo contains the
# correct license information in all files.
#
# script must be run from the $GIT_ROOT
#
# author: Mauricio Carneiro
# date: 1/9/13

echo "Checking all licenses in the GATK.. ";

result=0
ls public/**/*.java     | python private/python/licensing/CheckLicense.py
if [[ "$?" == "255" ]]
then
    result=1
fi

ls protected/**/*.java  | python private/python/licensing/CheckLicense.py
if [[ "$?" == "255" ]]
then
    result=1
fi

ls private/**/*.java    | python private/python/licensing/CheckLicense.py
if [[ "$?" == "255" ]]
then
    result=1
fi

ls public/**/*.scala    | python private/python/licensing/CheckLicense.py
if [[ "$?" == "255" ]]
then
    result=1
fi

ls protected/**/*.scala | python private/python/licensing/CheckLicense.py
if [[ "$?" == "255" ]]
then
    result=1
fi

ls private/**/*.scala   | python private/python/licensing/CheckLicense.py
if [[ "$?" == "255" ]]
then
    result=1
fi

exit $result

