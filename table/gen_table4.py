from preset import DATA_FILE_PATH
from preset import load_csv
from preset import dump_json
from preset import dump_csv
from collections import defaultdict

SHOW_VARIANT = [
    'B.1.1.7',
    'B.1.351',
    'P.1',
    'CAL.20C',
    'E484K',
    'N501Y',
    'K417N',
    'N439K',
    'Y453F',
    '∆69/70',
    '∆144',
]


def group_variants(records):
    variant_groups = defaultdict(list)
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
    cp_groups = defaultdict(list)
    for r in records:
        plasma = r['Plasma']
        if plasma.startswith('CP'):
            cp_groups['cp'].append(r)
        elif plasma.lower() in ('mild', 'severe'):
            cp_groups['cp'].append(r)
        else:
            cp_groups['vac'].append(r)

    result = {'variant': variant}
    for plasma, rec_list in cp_groups.items():
        num_studies = len(
            set([r['Reference'].replace('*', '') for r in rec_list]))
        num_s = sum([int(r['S']) for r in rec_list])
        num_i = sum([int(r['I']) for r in rec_list])
        num_r = sum([int(r['R']) for r in rec_list])
        num_samples = num_i + num_r + num_s

        pcnt_s = round(num_s / num_samples * 100)
        pcnt_i = round(num_i / num_samples * 100)
        pcnt_r = 100 - pcnt_s - pcnt_i

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


def gen_table4():
    cp_variant_records = load_csv(DATA_FILE_PATH / 'TableS4.csv')
    cp_mut_records = load_csv(DATA_FILE_PATH / 'TableS6.csv')

    variant_groups = {}
    variant_groups.update(group_variants(cp_variant_records))
    variant_groups.update(group_variants(cp_mut_records))

    result = []
    for variant in SHOW_VARIANT:
        records = variant_groups[variant]
        result.append(
            process_record(variant, records)
        )

    save_file = DATA_FILE_PATH / 'table4.json'
    dump_json(save_file, result)

    save_file = DATA_FILE_PATH / 'table4.csv'
    dump_csv(save_file, result)

    figure_results = []
    for item in result:
        variant = item['variant']
        for plasma in ['vac', 'cp']:
            for susc in ['s', 'i', 'r']:
                figure_results.append({
                    'variant': variant,
                    'study': plasma,
                    'num_study': item['{}_studies'.format(plasma)],
                    'susc': susc.upper(),
                    'sample_num': item['{}_num_{}_fold'.format(plasma, susc)]
                })

    save_file = DATA_FILE_PATH / 'table4-figure.csv'
    dump_csv(save_file, figure_results)
