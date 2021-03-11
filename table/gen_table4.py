from preset import DATA_FILE_PATH
from preset import load_csv
from preset import dump_json
from preset import dump_csv
from collections import defaultdict

SHOW_STRAIN = [
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


def group_strains(records):
    strain_groups = defaultdict(list)
    for r in records:
        strain = r['Strain name']
        if strain not in SHOW_STRAIN:
            continue
        strain_groups[strain].append(r)

    return strain_groups


def parse_fold(rec):
    fold = rec['Fold']
    if fold[0].isdigit():
        return float(fold)
    else:
        return float(fold[1:])


def process_record(strain, records):
    cp_groups = defaultdict(list)
    for r in records:
        plasma = r['Plasma']
        if plasma.startswith('CP'):
            cp_groups['cp'].append(r)
        else:
            cp_groups['vac'].append(r)

    result = {'strain': strain}
    for plasma, rec_list in cp_groups.items():
        num_studies = len(set([r['Reference'] for r in rec_list]))
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
    cp_strain_records = load_csv(DATA_FILE_PATH / 'TableS4.csv')
    cp_mut_records = load_csv(DATA_FILE_PATH / 'TableS6.csv')

    strain_groups = {}
    strain_groups.update(group_strains(cp_strain_records))
    strain_groups.update(group_strains(cp_mut_records))

    result = []
    for strain in SHOW_STRAIN:
        records = strain_groups[strain]
        result.append(
            process_record(strain, records)
        )

    save_file = DATA_FILE_PATH / 'table4.json'
    dump_json(save_file, result)

    save_file = DATA_FILE_PATH / 'table4.csv'
    dump_csv(save_file, result)
