#!/bin/zsh

# This script will update all files in the GATK with the license information
# in the files living in the license directory. This script is meant to be
# run manually (and hopefully only once). To run it, you must be in the
# $GATK root directory. If you want to run it from a different directory
# you will have to update the relative paths on this file
#
# author: Mauricio Carneiro
# date: 1/9/13

ls -1 public/**/*.java     | python private/python/licensing/UpdateLicense.py
ls -1 protected/**/*.java  | python private/python/licensing/UpdateLicense.py
ls -1 private/**/*.java    | python private/python/licensing/UpdateLicense.py
ls -1 public/**/*.scala    | python private/python/licensing/UpdateLicense.py
ls -1 protected/**/*.scala | python private/python/licensing/UpdateLicense.py
ls -1 private/**/*.scala   | python private/python/licensing/UpdateLicense.py
