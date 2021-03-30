from preset import DATA_FILE_PATH
from preset import load_csv
from preset import dump_json
from collections import defaultdict
from statistics import median
from preset import round_number

SHOW_VARIANT = [
    'B.1.1.7',
    'B.1.351',
    'P.1',
    'B.1.526',
    'B.1.429/B.1.427',
    'N501Y',
    'E484K',
    'K417N',
    'L452R',
    'N439K',
    'Y453F',
    'F490S',
    'S494P',
    # '∆69/70',
    # '∆144',
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
    'Casirivimab+Imdevimab': 'cas_imd',
    'Cilgavimab+Tixagevimab': 'cil_tix',
    'Bamlanivimab+Etesevimab': 'bam_ete',
    'BRII-196': 'BRII-196',
    'BRII-198': 'BRII-198',
    'C135': 'C135',
    'C144': 'C144',
    'JMB2002': 'JMB2002',
    'ADG20': 'ADG20',
}

DATA_PROBLEM = [
    ('B.1.1.7', 'sot'),
    ('N501Y', 'sot'),
]


def resistant_tester(x):
    return x >= 10


def partially_resistant_tester(x):
    return x > 3 and x < 10


COLOR_SETTING = {
    'resistant': {
        'tester': resistant_tester,
        'color': '#146aa8',
    },
    'partially-resistant': {
        'tester': partially_resistant_tester,
        'color': '#7fcbee'
    },
}


def get_color(medium_value):
    medium_value = float(medium_value)

    color = None
    for level, attr in COLOR_SETTING.items():
        tester = attr['tester']
        if tester(medium_value):
            color = attr['color']

    return color


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
        reference = reference.replace('*', '').replace('†', '')

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

        # rec_list = unique_reference(rec_list)

        rec_list.sort(key=parse_fold)
        fold_values = [
            100 if (i['Fold'] == '>100') else float(i['Fold'])
            for i in rec_list]
        medium_value = median(fold_values)

        if medium_value >= 100:
            medium_value_str = '&gt;100'
        else:
            medium_value_str = str(round_number(medium_value))
        # num_rec_list = len(set([row['Reference'] for i in rec_list]))
        num_rec_list = len(rec_list)

        tmpl = '{}<sub>{}</sub>'
        for s, m in DATA_PROBLEM:
            if s == variant and m == short_name:
                tmpl += '*'
                break

        result_markdown = tmpl.format(
            medium_value_str, num_rec_list if num_rec_list > 1 else ''
        )

        color = get_color(medium_value)
        if color:
            color_tmpl = ('<span '
                'style="padding: 0.1rem 0.2rem;margin: 0px 0.2rem;background-color:{color};color:white">'
                '{content}</span>')
            result_markdown = color_tmpl.format(
                color=color,
                content=result_markdown
            )

        result[short_name] = result_markdown

    return result


def gen_table_mab():
    mab_variant_records = load_csv(DATA_FILE_PATH / 'table_mab_variant.csv')
    mab_mut_records = load_csv(DATA_FILE_PATH / 'table_mab_muts.csv')

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

    save_file = DATA_FILE_PATH / 'table_mab.json'
    dump_json(save_file, result)
