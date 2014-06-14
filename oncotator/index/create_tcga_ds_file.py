# LICENSE_GOES_HERE
import sys
import pandas
import numpy


def add_full_mut_key_columns(input_df, chr_col_name='Chromosome', start_col_name='Start_position',
                             end_col_name='End_position', ref_col_name='Reference_Allele',
                             alt_col_name='Tumor_Seq_Allele1'):
    chrom_series, start_series = input_df[chr_col_name].apply(str), input_df[start_col_name].apply(str)
    end_series, ref_series = input_df[end_col_name].apply(str), input_df[ref_col_name]
    tum_series = input_df[alt_col_name]
    input_df[
        'full_mut_key'] = chrom_series + '_' + start_series + '_' + end_series + '_' + ref_series + '_' + tum_series


def main(tcga_mutations_fname, output_fname):
    print "Loading tcga data..."
    tcga_df = pandas.read_table(tcga_mutations_fname, skiprows=1, skipfooter=1)
    get_alt_allele = lambda x: x['tumor_seq_allele_1'] if x['tumor_seq_allele_1'] != x['reference_allele'] else x[
        'tumor_seq_allele_2']
    tcga_df['alt_allele'] = tcga_df.apply(get_alt_allele, axis=1)

    add_full_mut_key_columns(tcga_df, chr_col_name='chromosome', start_col_name='start_position',
                             end_col_name='end_position',
                             ref_col_name='reference_allele', alt_col_name='alt_allele')

    tcga_mut_counts = dict(tcga_df.full_mut_key.value_counts())
    all_tcga_diseases = set(tcga_df.tcga_disease)
    total_disease_sample_counts = dict()
    for tcga_disease in all_tcga_diseases:
        total_disease_sample_counts[tcga_disease] = len(
            set(tcga_df[tcga_df.tcga_disease == tcga_disease].tumor_sample_barcode))

    print "Indexing tcga data..."
    indexed_tcga_data = dict()
    for idx, row in tcga_df.iterrows():
        disease = row['tcga_disease']
        key = row['full_mut_key']
        if key not in indexed_tcga_data:
            indexed_tcga_data[key] = dict()

        if disease not in indexed_tcga_data[key]:
            indexed_tcga_data[key][disease] = 0

        indexed_tcga_data[key][disease] += 1

    print "Parsing tcga data..."
    tcga_mut_count = dict()
    tcga_muts_disease_names = dict()
    tcga_muts_disease_counts = dict()
    tcga_muts_disease_percents = dict()
    all_mut_keys = set()
    for mut_key in indexed_tcga_data:
        all_mut_keys.add(mut_key)
        disease_data = sorted(indexed_tcga_data[mut_key].items(), key=lambda x: x[1], reverse=True)
        tcga_muts_disease_names[mut_key] = ':'.join([d[0] for d in disease_data])
        tcga_muts_disease_counts[mut_key] = ':'.join([str(d[1]) for d in disease_data])
        get_mut_percent = lambda d: '%.2f' % (float(d[1]) / total_disease_sample_counts[d[0]] * 100)
        tcga_muts_disease_percents[mut_key] = ':'.join([get_mut_percent(d) for d in disease_data])
        tcga_mut_count[mut_key] = sum(d[1] for d in disease_data)

    print "Writing parsed tcga data..."
    output_headers = ['chr', 'start', 'end', 'ref_allele', 'alt_alele', 'tcga_mut_count', 'tcga_muts_disease_names',
                      'tcga_muts_disease_counts', 'tcga_muts_disease_percents']
    out_fh = open(output_fname, 'w')
    out_fh.write('\t'.join(output_headers) + '\n')

    all_mut_keys = sorted(all_mut_keys)
    for mut_key in all_mut_keys:
        output_fields = mut_key.split('_')
        output_fields.extend(
            [str(tcga_mut_count[mut_key]), tcga_muts_disease_names[mut_key], tcga_muts_disease_counts[mut_key],
             tcga_muts_disease_percents[mut_key]])
        out_fh.write('\t'.join(output_fields) + '\n')

    out_fh.close()


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "Usage: python create_tcga_ds_file.py input_tcga_file output_oncotator_ds_file"
        print "Example: create_tcga_ds_file.py mutations.2013042800.txt tcga_2013042800.tsv"
        sys.exit(1)

    tcga_mutations_fname, output_fname = sys.argv[1:]
    main(tcga_mutations_fname, output_fname)



