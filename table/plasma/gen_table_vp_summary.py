from preset import DATA_FILE_PATH
from preset import load_csv
from preset import dump_json
from preset import dump_csv

from .preset import IGNORE_VACCINE_NAME
from .common import get_num_fold_results_number_pair
from collections import defaultdict

SHOW_VARIANT = [
    'Alpha',
    'Beta',
    'Gamma',
    'Delta',
    'Omicron',
    'Iota',
    'Epsilon',
    'Kappa',
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
        variant = r['pattern']
        if variant not in SHOW_VARIANT:
            continue
        variant_groups[variant].append(r)

    return variant_groups


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
            set([r['ref_name'].replace('*', '')
                for r in rec_list]))

        aggre_list = [r for r in rec_list if 'Fold' in r.keys()]
        indiv_list = [r for r in rec_list if 'Fold' not in r.keys()]

        num_s, num_i, num_r = get_num_fold_results_number_pair(
            indiv_list, aggre_list)

        num_fold_results = num_i + num_r + num_s

        if num_fold_results:
            pcnt_s = round(num_s / num_fold_results * 100)
            pcnt_i = round(num_i / num_fold_results * 100)
            pcnt_r = 100 - pcnt_s - pcnt_i
        else:
            pcnt_i = 0
            pcnt_r = 0
            pcnt_s = 0

        result['vp_studies'] = num_studies
        result['vp_num_fold'] = num_fold_results
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
        if vaccine in IGNORE_VACCINE_NAME:
            continue
        vaccine_group[vaccine].append(rec)

    return vaccine_group


def gen_table_vp_summary():
    vp_variant_records = load_csv(DATA_FILE_PATH / 'table_vp_variants.csv')
    vp_mut_records = load_csv(DATA_FILE_PATH / 'table_vp_muts.csv')

    variant_groups = defaultdict(list)
    vp_variant_records = [
        v for v in vp_variant_records
        if (
            (int(float(v['dosage'])) != 1 and v['Plasma'] != 'Ad26.COV2.S')
            or
            (int(float(v['dosage'])) == 1 and v['Plasma'] == 'Ad26.COV2.S')
        )
    ]
    group_variants(variant_groups, vp_variant_records)
    vp_mut_records = [
        v for v in vp_mut_records
        if (
            (int(float(v['dosage'])) != 1 and v['Plasma'] != 'Ad26.COV2.S')
            or
            (int(float(v['dosage'])) == 1 and v['Plasma'] == 'Ad26.COV2.S')
        )
    ]
    group_variants(variant_groups, vp_mut_records)

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
