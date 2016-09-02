
**The docker scripts in here are meant for oncotator devs, not end users.  End users can get dockerhub images from https://hub.docker.com/r/broadinstitute/oncotator**

### I can't find the oncotator executable in the docker image.  Where is it?

By default, ``${ONCOTATOR}`` is ``/root/oncotator_venv/bin/oncotator``.  However, this is added to the path in the docker image, so ``oncotator ...`` should work.

### Where is the transcript override file in the docker image?

``/root/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt``

### How do I build a docker image based on a git hash?
``docker build -t broadinstitute/oncotator:${HASH} --build-arg ONCOTATOR_VERSION=${HASH} --build-arg GITHUB_DIR=\  .``

### What is Dockerfile.eval-gatk-protected ? 

This is provided as a courtesy to the GATK dev team for a very specific evaluation internal to the Broad.  ``example_input_targets.tsv`` is part of this image.

### Does the docker image include the default datasources?

No.  These must be added manually when running the docker image.  See docker docs for mounting volumes.

### Other notes:
- Oncotator is inside a venv on the docker image.
- Oncotator devs must still manually push images to dockerhub. 