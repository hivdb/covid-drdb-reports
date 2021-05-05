from preset import DATA_FILE_PATH
from preset import load_csv
from preset import dump_json
from preset import dump_csv

from .preset import VP_IGNORE
from .preset import VP_RENAME
from .common import get_sample_number_pair
from collections import defaultdict

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


def process_record(variant, records):

    vaccine_group = group_vaccine(records)

    result_list = {
        'variant': variant,
        'vaccine': []
    }

    for vaccine, rec_list in vaccine_group.items():

        result = {
            'name': vaccine
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

        result['vp_studies'] = num_studies
        result['vp_samples'] = num_samples
        result['vp_s_fold'] = '{}%'.format(
            pcnt_s) if pcnt_s else 0
        result['vp_i_fold'] = '{}%'.format(
            pcnt_i) if pcnt_i else 0
        result['vp_r_fold'] = '{}%'.format(
            pcnt_r) if pcnt_r else 0

        result['vp_num_s_fold'] = num_s
        result['vp_num_i_fold'] = num_i
        result['vp_num_r_fold'] = num_r

        result_list['vaccine'].append(result)

    return result_list


def group_vaccine(rx_list):

    vaccine_group = defaultdict(list)
    for rec in rx_list:
        vaccine = rec['Plasma']
        ref_name = rec['Reference']
        if (ref_name, vaccine) in VP_IGNORE:
            continue

        vaccine_name = VP_RENAME.get(vaccine)
        if not vaccine_name:
            print('Unknown vaccine', ref_name, vaccine)
            continue

        vaccine_group[vaccine].append(rec)

    return vaccine_group


def gen_table_vp_summary():
    cp_variant_records = load_csv(DATA_FILE_PATH / 'table_vp_variants.csv')
    cp_mut_records = load_csv(DATA_FILE_PATH / 'table_vp_muts.csv')

    variant_groups = defaultdict(list)
    group_variants(variant_groups, cp_variant_records)
    group_variants(variant_groups, cp_mut_records)

    results = []
    for variant in SHOW_VARIANT:
        records = variant_groups.get(variant, [])

        results.append(
            process_record(variant, records)
        )

    save_file = DATA_FILE_PATH / 'table_vp_summary.json'
    dump_json(save_file, results)

    csv_result = []
    for rec in results:
        variant = rec['variant']
        for item in rec['vaccine']:
            item['variant'] = variant
            csv_result.append(item)

    save_file = DATA_FILE_PATH / 'table_vp_summary.csv'
    dump_csv(save_file, csv_result)
