from preset import load_csv
from preset import dump_csv
from preset import dump_json
from preset import DATA_FILE_PATH
from scipy.stats.mstats import gmean
from preset import round_number
from preset import group_records_by


WT_LIST = [
    'A',
    'A.1',
    'A.2.2',
    'B.1',
    'B',
    'B.1.2',
    'B.1.1.117',
    'B.1.177',
    'B.1.319',
    'B.1.1',
]


VARIANT_LIST = [
    'WT',
    'Alpha',
    'Beta',
    'Gamma',
    'Delta',
    'Omicron/BA.1',
    'Omicron/BA.1.1',
    'Omicron/BA.2',
    'Omicron/BA.2.12.1',
    'Omicron/BA.3',
    'Omicron/BA.4/5',
    'Omicron/BA.2.75',
]


def gen_cross_neutral(
        src=DATA_FILE_PATH / 'figure' / 'figure_plasma_variant_titer.csv',
        dst=DATA_FILE_PATH / 'plasma' / 'cross_neutral.csv'):

    table = load_csv(src)

    table = [
        i
        for i in table
        if not i['dosage'] and i['infection']
    ]

    [
        i.update({
            'infection': 'WT' if i['infection'] in WT_LIST else i['infection'],
            'var_name': 'WT' if i['var_name'] in WT_LIST else i['var_name']
        })
        for i in table
    ]

    # print(set([
    #     i['var_name'] for i in table if i['var_name'] not in VARIANT_LIST]))

    table = [
        i
        for i in table
        if i['infection'] in VARIANT_LIST
        and i['var_name'] in VARIANT_LIST
    ]

    infection_variants = [
        i
        for i in VARIANT_LIST
        if i in [
            j['infection']
            for j in table
        ]
    ]
    results = []

    get_result_table(table, infection_variants, results)

    dump_csv(dst, results)

    dst = dst.with_suffix('.json').with_stem(f'table_{dst.stem}')
    dst = dst.parent.parent / dst.name
    results = get_json_results(results)
    dump_json(dst, results)


def get_result_table(table, infection_variants, results):
    for test in VARIANT_LIST:
        for inf in infection_variants:
            rec_list = [
                i
                for i in table
                if i['infection'] == inf and i['var_name'] == test
            ]

            get_result_by_month(test, inf, rec_list, results, 0, 1)
            get_result_by_month(test, inf, rec_list, results, 1, 6)
            get_result_by_month(test, inf, rec_list, results, 6, 1000)


def get_result_by_month(test, inf, rec_list, results, start_month, stop_month):
    rec_list = [
        i
        for i in rec_list
        if i['month'] > start_month and i['month'] <= stop_month
    ]

    titer_list = [i['titer'] for i in rec_list]
    num_exp_list = [i['num_exp'] for i in rec_list]
    num_ref_name = len(set(
        [i['ref_name'] for i in rec_list]
    ))
    if not titer_list or not num_exp_list:
        results.append({
            'infection': inf,
            'test': test,
            'geomean': '',
            'num_ref_name': num_ref_name,
            'group_name': f'{start_month}-{stop_month}'
        })
        return

    geomean = gmean(
        titer_list,
        weights=num_exp_list
    )
    results.append({
        'infection': inf,
        'test': test,
        'geomean': round_number(geomean),
        'num_ref_name': num_ref_name,
        'group_name': f'{start_month}-{stop_month}'
    })


def get_json_results(table):

    results = []

    for test, rec_list in group_records_by(table, 'test').items():
        rec = {
            'test': test.replace('Omicron/', '')
        }
        for inf, inf_rec_list in group_records_by(rec_list, 'infection').items():
            for idx, group_name in enumerate([
                    '0-1', '1-6', '6-1000']):
                group_data = [
                    i
                    for i in inf_rec_list
                    if i['group_name'] == group_name
                    ][0]

                rec[f'{inf}_{idx+1}'] = {
                    'geomean': group_data['geomean'],
                    'num_ref_name': group_data['num_ref_name']
                }

        results.append(rec)

    return results
