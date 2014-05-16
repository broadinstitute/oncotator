#!/bin/zsh

if [[ "$2" == "-debug" ]]
then
    echo $1
fi

if [[ ! ("${1##*.}" == "java" || "${1##*.}" == "scala") ]]
then
    echo "Only java or scala files should be updated for license. This is a $1:e file";
    exit 1;
fi

isPrivate=`echo "$1" | sed -n '/private/p'`;
isProtected=`echo "$1" | sed -n '/protected/p'`;
isPublic=`echo "$1" | sed -n '/public/p'`;

if [[ $isPrivate != "" ]]
then
    licenseFile=licensing/private_license.txt;
elif [[ $isProtected != "" ]]
then
    licenseFile=licensing/protected_license.txt;
elif [[ $isPublic != "" ]]
then
    licenseFile=licensing/public_license.txt;
else
    echo "$1 is not in public, private or protected and is a source file. Please place it in the appropriate directory so we can choose the appropriate license";
    exit 1;
fi

python private/python/UpdateLicense.py $1 $licenseFile > $1.newlicense && mv $1.newlicense $1;
