from preset import DATA_FILE_PATH
from preset import load_csv
from preset import dump_json
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
    'L452R',
    'Y453F',
    'F490S',
    'S494P',
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


# def get_mab_names(mab_name):
#     return [mab_name]
#     if '/' in mab_name:
#         mab_names = mab_name.split('/')
#     elif '+' in mab_name:
#         mab_names = mab_name.split('+')
#     else:
#         mab_names = mab_name.split()
#     mab_names = [m.strip() for m in mab_names]

#     return mab_names


def unique_reference(rec_list):
    unique_rec_list = {}
    for r in rec_list:
        reference = r['Reference']
        reference = reference.replace('*', '')

        previous = unique_rec_list.get(reference)
        if not previous:
            unique_rec_list[reference] = r
            continue

        if previous['Fold'] < r['Fold']:
            unique_rec_list[reference] = r

    rec_list = list(unique_rec_list.values())

    return rec_list


def process_record(variant, records):
    mab_groups = defaultdict(list)
    for r in records:
        mab_name = r['Mab name']
        if mab_name in SHOW_MABS.keys():
            short_name = SHOW_MABS[mab_name]
            mab_groups[short_name].append(r)

    result = {'variant': variant}
    for mab_name, short_name in SHOW_MABS.items():
        rec_list = mab_groups.get(short_name)

        if not rec_list:
            result[short_name] = '-'
            continue

        rec_list = unique_reference(rec_list)

        rec_list.sort(key=parse_fold)
        max_value = rec_list[-1]['Fold']
        max_value = max_value.replace('>', '&gt;')
        num_rec_list = len(set([i['Reference'] for i in rec_list]))

        tmpl = '{}<sub>{}</sub>'
        for s, m in DATA_PROBLEM:
            if s == variant and m == short_name:
                tmpl += '*'
                break

        result[short_name] = tmpl.format(
            max_value, num_rec_list if num_rec_list > 1 else ''
        )

    return result


def gen_table3():
    mab_variant_records = load_csv(DATA_FILE_PATH / 'TableS5.csv')
    mab_mut_records = load_csv(DATA_FILE_PATH / 'TableS7.csv')

    variant_groups = {}
    variant_groups.update(group_variants(mab_variant_records))
    variant_groups.update(group_variants(mab_mut_records))

    result = []
    for variant in SHOW_VARIANT:
        records = variant_groups.get(variant)
        if not records:
            continue
        result.append(
            process_record(variant, records)
        )

    save_file = DATA_FILE_PATH / 'table3.json'
    dump_json(save_file, result)
