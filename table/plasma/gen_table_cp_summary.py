from preset import DATA_FILE_PATH
from preset import load_csv
from preset import dump_json
from preset import dump_csv

from collections import defaultdict
from .common import get_sample_number_pair

SHOW_VARIANT = [
    'B.1.1.7',
    'B.1.351',
    'P.1',
    'B.1.526',
    'B.1.427/9',
    'B.1.617',
    'N501Y',
    'E484K',
    'K417N',
    'L452R',
    'N439K',
    'Y453F',
    # 'F490S',
    # 'S494P',
    '∆69/70',
    '∆144',
]


def group_variants(variant_groups, records):
    for r in records:
        variant = r['Variant name']
        if variant not in SHOW_VARIANT:
            continue
        variant_groups[variant].append(r)

    return variant_groups


def parse_fold(rec):
    fold = rec['Fold']
    if fold[0].isdigit():
        return float(fold)
    else:
        return float(fold[1:])


def process_record(variant, rec_list):

    result = {
        'variant': variant
    }

    num_studies = len(
        set([r['Reference'].replace('*', '')
            for r in rec_list]))

    aggre_list = [r for r in rec_list if 'Fold' in r.keys()]
    indiv_list = [r for r in rec_list if 'Fold' not in r.keys()]

    num_s, num_i, num_r = get_sample_number_pair(indiv_list, aggre_list)

    num_samples = num_i + num_r + num_s

    if num_samples:
        pcnt_s = round(num_s / num_samples * 100)
        pcnt_i = round(num_i / num_samples * 100)
        pcnt_r = 100 - pcnt_s - pcnt_i
    else:
        pcnt_i = 0
        pcnt_r = 0
        pcnt_s = 0

    result['cp_studies'] = num_studies
    result['cp_samples'] = num_samples
    result['cp_s_fold'] = '{}%'.format(
        pcnt_s) if pcnt_s else 0
    result['cp_i_fold'] = '{}%'.format(
        pcnt_i) if pcnt_i else 0
    result['cp_r_fold'] = '{}%'.format(
        pcnt_r) if pcnt_r else 0

    result['cp_num_s_fold'] = num_s
    result['cp_num_i_fold'] = num_i
    result['cp_num_r_fold'] = num_r

    return result


def gen_table_cp_summary():
    cp_variant_records = load_csv(DATA_FILE_PATH / 'table_cp_variants.csv')
    cp_mut_records = load_csv(DATA_FILE_PATH / 'table_cp_muts.csv')

    variant_groups = defaultdict(list)
    group_variants(variant_groups, cp_variant_records)
    group_variants(variant_groups, cp_mut_records)

    result = []
    for variant in SHOW_VARIANT:
        records = variant_groups.get(variant, [])
        result.append(
            process_record(variant, records)
        )

    save_file = DATA_FILE_PATH / 'table_cp_summary.json'
    dump_json(save_file, result)

    save_file = DATA_FILE_PATH / 'table_cp_summary.csv'
    dump_csv(save_file, result)
