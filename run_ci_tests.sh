#!/bin/bash
# Have script stop if there is an error
set -e

echo "This script must be run from the same directory as setup.py"

VENV=oncotator_nosetest_env_auto
bash ./scripts/create_oncotator_venv.sh $VENV
source $VENV/bin/activate
nosetests --all-modules --exe --with-xunit -w test -v --processes=4 --process-timeout=480  --process-restartworker
deactivate
rm -Rf $VENV