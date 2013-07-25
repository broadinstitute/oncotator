#!/bin/bash
# Have script stop if there is an error
set -e

echo "This script must be run from the same directory as setup.py"

VENV=oncotator_nosetest_env_auto
bash ./scripts/create_oncotator_venv.sh $VENV
source $VENV/bin/activate

# nosetests will return a non-zero error code if any unit test fails, so we do not want to stop this script
#   So we remove the error catching temporarily
set +e
nosetests --all-modules --exe --with-xunit -w test -v --processes=4 --process-timeout=480  --process-restartworker
set -e

deactivate
rm -Rf $VENV