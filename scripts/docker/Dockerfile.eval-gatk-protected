FROM broadinstitute/oncotator:1.9.0.0

RUN wget www.broadinstitute.org/~lichtens/ref_hg_gencode.tar.gz
RUN tar zxvf ref_hg_gencode.tar.gz
RUN rm -f ref_hg_gencode.tar.gz
RUN ln -s /root/xchip/cga/reference/annotation/db/ref_hg_gencode/ /root/onco_dbdir
COPY example_input_targets.tsv /root/eval-gatk-protected/example_input_targets.tsv

# Create a cache for CRSP targets that should speed the oncotator results
RUN /root/oncotator_venv/bin/oncotator --db-dir /root/onco_dbdir/ -c /root/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt -u file:///root/onco_cache/ -v /root/eval-gatk-protected/example_input_targets.tsv /root/eval-gatk-protected/TEST.per_target.oncotated.txt hg19 -i SEG_FILE -o SIMPLE_TSV
