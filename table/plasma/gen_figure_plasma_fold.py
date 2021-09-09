from preset import DATA_FILE_PATH
from preset import load_csv
from preset import dump_csv
from .preset import IGNORE_VACCINE_NAME

from .common import get_sample_number_pair
from collections import defaultdict

SHOW_VARIANT = [
    'Alpha',
    'Beta',
    'Gamma',
    'Delta',
    'Iota',
    'Epsilon',
    'Kappa',
    'N501Y',
    'E484K',
    'L452R',
    # 'K417N',
    # 'N439K',
    # 'Y453F',
    # 'F490S',
    # 'S494P',
    # '∆69/70',
    # '∆144',
]


REGULAR_PLASMA_NAMES = ['CP', 'BNT162b2', 'mRNA-1273', 'AZD1222']


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
    global REGULAR_PLASMA_NAMES

    cp_groups = defaultdict(list)

    for r in records:
        plasma = r['Plasma']
        if plasma.startswith('CP'):
            cp_groups['CP'].append(r)
        elif plasma.lower() in ('mild', 'severe'):
            cp_groups['CP'].append(r)
        elif plasma == 'IVIG':
            cp_groups['CP'].append(r)
        elif plasma in IGNORE_VACCINE_NAME:
            continue
        else:
            if plasma not in REGULAR_PLASMA_NAMES:
                REGULAR_PLASMA_NAMES.append(plasma)
            cp_groups[plasma].append(r)

    for name in REGULAR_PLASMA_NAMES:
        if name not in cp_groups:
            cp_groups['vac'] = []

    result = {'variant': variant}
    for plasma, rec_list in cp_groups.items():
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

        result['{}_studies'.format(plasma)] = num_studies
        result['{}_samples'.format(plasma)] = num_samples
        result['{}_s_fold'.format(plasma)] = '{}%'.format(
            pcnt_s) if pcnt_s else 0
        result['{}_i_fold'.format(plasma)] = '{}%'.format(
            pcnt_i) if pcnt_i else 0
        result['{}_r_fold'.format(plasma)] = '{}%'.format(
            pcnt_r) if pcnt_r else 0

        result['{}_num_s_fold'.format(plasma)] = num_s
        result['{}_num_i_fold'.format(plasma)] = num_i
        result['{}_num_r_fold'.format(plasma)] = num_r

    return result


def gen_figure_plasma_fold():
    cp_variant_records = load_csv(DATA_FILE_PATH / 'table_cp_variants.csv')
    cp_mut_records = load_csv(DATA_FILE_PATH / 'table_cp_muts.csv')
    vp_variant_records = load_csv(DATA_FILE_PATH / 'table_vp_variants.csv')
    vp_mut_records = load_csv(DATA_FILE_PATH / 'table_vp_muts.csv')

    variant_groups = defaultdict(list)
    group_variants(variant_groups, cp_variant_records)
    group_variants(variant_groups, cp_mut_records)

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

    result = []
    for variant in SHOW_VARIANT:
        records = variant_groups.get(variant, [])
        result.append(
            process_record(variant, records)
        )

    figure_results = []
    for item in result:
        variant = item['variant']
        for plasma in REGULAR_PLASMA_NAMES:
            for susc in ['s', 'i', 'r']:
                figure_results.append({
                    'variant': variant,
                    'study': plasma,
                    'num_study': item.get('{}_studies'.format(plasma), 0),
                    'susc': susc.upper(),
                    'sample_num': item.get(
                        '{}_num_{}_fold'.format(plasma, susc), 0)
                })

    save_file = DATA_FILE_PATH / 'figure_plasma_fold.csv'
    dump_csv(save_file, figure_results)
