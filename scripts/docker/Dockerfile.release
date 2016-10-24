FROM python:2.7.9

MAINTAINER oncotator

# This can be overridden with ``--build-arg ONCOTATOR_VERSION=<desired_version_here>``
ARG ONCOTATOR_VERSION=1.9.1.0

# If you are overriding the version and wish to build against a git hash, then add ``--build-arg GITHUB_DIR=\ ``
ARG GITHUB_DIR=tags/

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y wget curl unzip python-pip libxft-dev libfreetype6 libfreetype6-dev gcc git

#### Specific for google cloud support
RUN wget https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.zip \
    && unzip google-cloud-sdk.zip \
    && rm google-cloud-sdk.zip

RUN google-cloud-sdk/install.sh --usage-reporting=true --path-update=true --bash-completion=true --rc-path=/.bashrc --disable-installation-options
VOLUME ["/root/.config"]
ENV PATH /google-cloud-sdk/bin:$PATH

RUN yes | gcloud components update
RUN yes | gcloud components update preview
###########

## Add more python packages
RUN apt-get install -y python-dev python-setuptools emacs less lynx
RUN pip install virtualenv
###################


ENV ONCOTATOR_VERSION ${ONCOTATOR_VERSION}
ENV ONCOTATOR_VENV /root/oncotator_venv/
ENV GITHUB_DIR ${GITHUB_DIR}

WORKDIR /root/
RUN git clone https://github.com/broadinstitute/oncotator.git

WORKDIR /root/oncotator/
RUN git checkout ${GITHUB_DIR}v${ONCOTATOR_VERSION}
RUN bash scripts/create_oncotator_venv.sh -e ${ONCOTATOR_VENV}
RUN ${ONCOTATOR_VENV}/bin/python setup.py install
ENV ONCOTATOR ${ONCOTATOR_VENV}/bin/oncotator
ENV PATH ${ONCOTATOR_VENV}/bin/oncotator:$PATH

RUN ${ONCOTATOR} --help

# Get transcript override file
WORKDIR /root/
COPY tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt .

CMD ["bash"]

# END