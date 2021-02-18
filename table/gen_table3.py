from preset import DATA_FILE_PATH
from preset import load_csv
from preset import dump_json
from collections import defaultdict

SHOW_STRAIN = [
    'B.1.1.7',
    'B.1.351',
    'P.1',
    'E484K',
    'N501Y',
    'K417N',
    'N439K',
    'L452R',
    'Y453F',
]

SHOW_MABS = {
    'Casirivimab': 'cas',
    'Etesevimab': 'ete',
    'Tixagevimab': 'tix',
    'Bamlanivimab': 'bam',
    'Cilgavimab': 'cil',
    'Imdevimab': 'imd',
    'Sotrovimab': 'sot',
    'Regdanvimab': 'reg',
    'BRII-196': 'BRII-196',
    'BRII-198': 'BRII-198',
    'C135': 'C135',
    'C144': 'C144',
    'JMB2002': 'JMB2002',
    'S2E12': 'S2E12',
}

DATA_PROBLEM = [
    ('B.1.1.7', 'sot'),
    ('N501Y', 'sot'),
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
    mab_groups = defaultdict(list)
    for r in records:
        mab_name = r['Mab name']
        if mab_name not in SHOW_MABS.keys():
            continue
        mab_groups[mab_name].append(r)

    result = {'strain': strain}
    for mab_name, short_name in SHOW_MABS.items():
        rec_list = mab_groups.get(mab_name)
        if not rec_list:
            result[short_name] = '-'
            continue

        rec_list.sort(key=parse_fold)
        max_value = rec_list[-1]['Fold']
        max_value = max_value.replace('>', '&gt;')
        num_rec_list = len(set([i['Reference'] for i in rec_list]))

        tmpl = '{}<sub>{}</sub>'
        for s, m in DATA_PROBLEM:
            if s == strain and m == short_name:
                tmpl += '*'
                break

        result[short_name] = tmpl.format(
            max_value, num_rec_list if num_rec_list > 1 else ''
        )

    return result


def gen_table3():
    mab_strain_records = load_csv(DATA_FILE_PATH / 'TableS5.csv')
    mab_mut_records = load_csv(DATA_FILE_PATH / 'TableS7.csv')

    strain_groups = {}
    strain_groups.update(group_strains(mab_strain_records))
    strain_groups.update(group_strains(mab_mut_records))

    result = []
    for strain in SHOW_STRAIN:
        records = strain_groups[strain]
        result.append(
            process_record(strain, records)
        )

    save_file = DATA_FILE_PATH / 'table3.json'
    dump_json(save_file, result)
