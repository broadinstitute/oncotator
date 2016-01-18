FROM python:2

MAINTAINER Alex Ramos <aramos@broadinstitute.org>

RUN apt-get update && apt-get install unzip

RUN pip install numpy

RUN wget --no-check-certificate 'https://github.com/elephanthunter/PyVCF/archive/master.zip'

RUN unzip master.zip && cd PyVCF-master && python setup.py install && cd .. && rm -Rf PyVCF-master && rm -f master.*

RUN pip install ngslib

ADD . /oncotator

RUN cd oncotator/ && \
wget https://pysam.googlecode.com/files/pysam-0.7.5.tar.gz && \
python /oncotator/distribute_setup.py && \
easy_install -U distribute && \
pip install pysam-0.7.5.tar.gz && \
python setup.py install

ENTRYPOINT ["Oncotator"]

CMD ["-h"]

# EXAMPLE BUILD CMD
# docker build -t oncotator .

# EXAMPLE RUN CMDS
# docker run -it oncotator -h
# docker run -v /path/to/data:/data -v /path/to/oncotator_db_dir:/db_dir -it oncotator --db-dir /db_dir /data/in.maflite /data/out.maf.txt hg19
