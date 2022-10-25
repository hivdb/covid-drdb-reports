from pathlib import Path
from preset import DATA_FILE_PATH
from preset import load_csv
from preset import load_yaml
from preset import dump_json
from collections import defaultdict
from statistics import median
from resistancy import round_fold
from resistancy import parse_fold
from z_score import get_outlier

from resistancy import is_partial_resistant
from resistancy import is_resistant

CONFIG = load_yaml(Path(__file__).resolve().parent / 'mab_report.yml')

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


def group_variants(records, show_list):
    variant_groups = defaultdict(list)
    for r in records:
        variant = r['pattern']
        if variant not in show_list:
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

    result = {'variant': variant}

    if not records:
        for mab_name, short_name in CONFIG['SHOW_MABS'].items():
            result[short_name] = {}
            result[short_name]['fold'] = '-'

        return result

    for r in records:
        mab_name = r['ab_name']
        if mab_name in CONFIG['SHOW_MABS'].keys():
            short_name = CONFIG['SHOW_MABS'][mab_name]
            mab_groups[short_name].append(r)

    for mab_name, short_name in CONFIG['SHOW_MABS'].items():
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

        # outlier = get_outlier(rec_list, 'fold')
        # if outlier:
        #     print('Outliers', outlier)

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


def mab_report(
        load_folder=DATA_FILE_PATH / 'mab',
        dump_folder=DATA_FILE_PATH):

    _mab_report(
        mab_variants_file=load_folder / 'table_mab_variant.csv',
        mab_muts_file=load_folder / 'table_mab_muts.csv',
        dump_file=dump_folder/'table_mab.json',
        show_list=CONFIG['WT_VARIANTS']
    )

    _mab_report(
        mab_variants_file=load_folder / 'table_mab_omicron_variant.csv',
        mab_muts_file=load_folder / 'table_mab_omicron_muts.csv',
        dump_file=dump_folder/'table_mab_ba2.json',
        show_list=CONFIG['OMICRON_VARIANTS']
    )


def _mab_report(mab_variants_file, mab_muts_file, dump_file, show_list):

    mab_variant_records = (
        load_csv(mab_variants_file) if mab_variants_file.exists() else []
    )
    mab_mut_records = (
        load_csv(mab_muts_file) if mab_muts_file.exists() else []
    )

    variant_groups = {}
    variant_groups.update(group_variants(mab_variant_records, show_list))
    variant_groups.update(group_variants(mab_mut_records, show_list))

    result = []
    for variant in show_list:
        records = variant_groups.get(variant, [])
        if not records and not variant.startswith('Omicron'):
            continue

        result.append(
            process_record(variant, records)
        )

    dump_json(dump_file, result)
