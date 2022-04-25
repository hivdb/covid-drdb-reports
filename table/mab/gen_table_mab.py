from preset import DATA_FILE_PATH
from preset import load_csv
from preset import dump_json
from collections import defaultdict
from statistics import median
from resistancy import round_fold
from resistancy import parse_fold
from z_score import get_outlier

from resistancy import is_partial_resistant
from resistancy import is_resistant

SHOW_VARIANT = [
    'Alpha',
    'Beta',
    'Gamma',
    'Delta',
    'Omicron/BA.1',
    'Omicron/BA.2',
    'Iota',
    'Epsilon',
    'Kappa',
    'Lambda',
    'Mu',
    'N501Y',
    'E484K',
    'K417N',
    'L452R',
    'T478K',
    'N439K',
    'Y453F',
    'F490S',
    'S494P',
    'S371F',
    'S371L',
    'D405N',
    'N440K',
    'G446S',
    'E484A',
    'Q493R',
    ''
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
    'Bebtelovimab': 'beb',
    'Casirivimab/Imdevimab': 'cas_imd',
    'Cilgavimab/Tixagevimab': 'cil_tix',
    'Bamlanivimab/Etesevimab': 'bam_ete',
    'C135': 'C135',
    'C144': 'C144',
    'C135/C144': 'c_135_144',
    'Amubarvimab': 'b_196',
    'Romlusevimab': 'b_198',
    'Amubarvimab/Romlusevimab': 'b_196_198',
    'JMB2002': 'JMB2002',
    'Adintrevimab': 'Adintrevimab',
    'Vir-7832': 'Vir7832',
    'DXP-604': 'DXP604',
}

# DATA_PROBLEM = [
#     ('Alpha', 'sot'),
#     ('N501Y', 'sot'),
# ]


COLOR_SETTING = {
    'resistant': {
        'tester': is_resistant,
        'color': '#146aa8',
    },
    'partially-resistant': {
        'tester': is_partial_resistant,
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
        variant = r['pattern']
        if variant not in SHOW_VARIANT:
            continue
        variant_groups[variant].append(r)

    return variant_groups


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
        ref_name = r['ref_name']
        ref_name = ref_name.replace('*', '').replace('†', '')

        previous = unique_rec_list.get(ref_name)
        if not previous:
            unique_rec_list[ref_name] = r
            continue

        if previous['fold'] < r['fold']:
            unique_rec_list[ref_name] = r

    rec_list = list(unique_rec_list.values())

    return rec_list


def process_record(variant, records):
    mab_groups = defaultdict(list)
    for r in records:
        mab_name = r['ab_name']
        if mab_name in SHOW_MABS.keys():
            short_name = SHOW_MABS[mab_name]
            mab_groups[short_name].append(r)

    result = {'variant': variant}
    for mab_name, short_name in SHOW_MABS.items():
        result[short_name] = {}
        rec_list = mab_groups.get(short_name)

        if not rec_list:
            result[short_name]['fold'] = '-'
            continue

        # rec_list = unique_reference(rec_list)

        [
            rec.update({
                'fold': parse_fold(rec['fold'])
            })
            for rec in rec_list
        ]

        rec_list.sort(key=lambda i: i['fold'])
        fold_values = [rec['fold'] for rec in rec_list]
        medium_value = median(fold_values)

        outlier = get_outlier(rec_list, 'fold')
        if outlier:
            print('Outliers', outlier)

        if medium_value >= 1000:
            fold_cmp = '>'
            fold = '1000'
            # medium_value_str = '&gt;1000'
        else:
            # medium_value_str = str(round_fold(medium_value))
            fold_cmp = ''
            fold = str(round_fold(medium_value))
        # num_rec_list = len(set([row['ref_name'] for i in rec_list]))
        num_rec_list = len(rec_list)

        # tmpl = '{}<sub>{}</sub>'
        # for s, m in DATA_PROBLEM:
        #     if s == variant and m == short_name:
        #         # tmpl += '*'
        #         fold = fold + '*'
        #         break

        result[short_name]['fold'] = fold
        result[short_name]['fold_cmp'] = fold_cmp
        result[short_name]['num'] = num_rec_list

        # result_markdown = tmpl.format(
        #     medium_value_str, num_rec_list if num_rec_list > 1 else ''
        # )

        # color = get_color(medium_value)
        # if color:
        #     color_tmpl = (
        #         '<span '
        #         'style="padding: 0.1rem 0.2rem;margin: '
        #         '0px 0.2rem;background-color:{color};color:white">'
        #         '{content}</span>')
        #     result_markdown = color_tmpl.format(
        #         color=color,
        #         content=result_markdown
        #     )

        # result[short_name] = result_markdown

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
