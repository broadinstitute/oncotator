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


def add_full_mut_key_to_cosmic(input_df):
    input_df.chromosome = input_df.chromosome.apply(lambda x: {'23': 'X', '24': 'Y', '25': 'MT'}.get(str(x), str(x)))
    chrom_series, start_series = input_df.chromosome, input_df['GRCh37 start'].apply(str)
    end_series = input_df['GRCh37 stop'].apply(str)
    ref_series = input_df.mut_nt.apply(lambda x: str(x).split('>')[0].upper())
    tum_series = input_df.mut_nt.apply(lambda x: str(x).split('>')[1].upper())
    input_df[
        'full_mut_key'] = chrom_series + '_' + start_series + '_' + end_series + '_' + ref_series + '_' + tum_series


def main(cosmic_mutations_fname, output_fname):
    print "Loading cosmic data..."
    cosmic_df = pandas.read_csv(cosmic_mutations_fname, skipinitialspace=True, skipfooter=3)
    cosmic_df = cosmic_df[cosmic_df.mut_nt.apply(lambda x: numpy.nan is not x)]
    add_full_mut_key_to_cosmic(cosmic_df)

    print "Indexing cosmic data..."
    indexed_cosmic_data = dict()
    for idx, row in cosmic_df.iterrows():
        disease = row['tumour_site']
        key = row['full_mut_key']
        if key not in indexed_cosmic_data:
            indexed_cosmic_data[key] = dict()

        if disease not in indexed_cosmic_data[key]:
            indexed_cosmic_data[key][disease] = [0, 0]

        indexed_cosmic_data[key][disease][0] += row['mutated_samples']
        indexed_cosmic_data[key][disease][1] += row['examined_samples']

    print "Parsing cosmic data..."
    cosmic_mut_count = dict()
    cosmic_muts_disease_names = dict()
    cosmic_muts_disease_counts = dict()
    cosmic_muts_disease_percents = dict()
    all_mut_keys = set()
    for mut_key in indexed_cosmic_data:
        all_mut_keys.add(mut_key)
        disease_data = sorted([(k, v[0], v[1]) for k, v in indexed_cosmic_data[mut_key].items()], key=lambda x: x[1],
                              reverse=True)
        cosmic_muts_disease_names[mut_key] = ':'.join([d[0] for d in disease_data])
        cosmic_muts_disease_counts[mut_key] = ':'.join([str(d[1]) for d in disease_data])
        get_mut_percent = lambda d: '%.2f' % (float(d[1]) / d[2] * 100)
        cosmic_muts_disease_percents[mut_key] = ':'.join([get_mut_percent(d) for d in disease_data])
        cosmic_mut_count[mut_key] = sum(d[1] for d in disease_data)

    print "Writing parsed cosmic data..."
    output_headers = ['chr', 'start', 'end', 'ref_allele', 'alt_alele', 'cosmic_mut_count', 'cosmic_muts_disease_names',
                      'cosmic_muts_disease_counts', 'cosmic_muts_disease_percents']
    out_fh = open(output_fname, 'w')
    out_fh.write('\t'.join(output_headers) + '\n')

    all_mut_keys = sorted(all_mut_keys)
    for mut_key in all_mut_keys:
        output_fields = mut_key.split('_')
        output_fields.extend(
            [str(cosmic_mut_count[mut_key]), cosmic_muts_disease_names[mut_key], cosmic_muts_disease_counts[mut_key],
             cosmic_muts_disease_percents[mut_key]])
        out_fh.write('\t'.join(output_fields) + '\n')

    out_fh.close()


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "Usage: python create_cosmic_ds_file.py input_cosmic_file output_oncotator_ds_file"
        print "Example: create_cosmic_ds_file.py UCSCMutExp_v65_100613.csv cosmic_v65.tsv"
        sys.exit(1)

    cosmic_mutations_fname, output_fname = sys.argv[1:]
    main(cosmic_mutations_fname, output_fname)



