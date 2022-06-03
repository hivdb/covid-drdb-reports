from unittest import result
from preset import load_csv
from preset import dump_csv
from preset import dump_json
from preset import DATA_FILE_PATH
from scipy.stats.mstats import gmean
from preset import round_number


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
    'Omicron/BA.4',
    'Omicron/BA.5',
]


def gen_cross_neutral(
        src=DATA_FILE_PATH / 'figure' / 'figure_vp_variant_titer.csv',
        dst=DATA_FILE_PATH / 'plasma' / 'cross_neutraling.csv'):

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
        if i['infection'] in VARIANT_LIST and i['var_name'] in VARIANT_LIST
    ]

    infection_variants = [
        i
        for i in VARIANT_LIST
        if i in [
            j['infection']
            for j in table
        ]
    ]

    json_results = []
    csv_results = []
    for test in VARIANT_LIST:
        rec = {
            'test': test
        }
        for inf in infection_variants:
            rec_list = [
                i
                for i in table
                if i['infection'] == inf and i['var_name'] == test
            ]
            titer_list = [i['titer'] for i in rec_list]
            num_exp_list = [i['num_exp'] for i in rec_list]
            num_ref_name = len(set(
                [i['ref_name'] for i in rec_list]
            ))
            if not titer_list or not num_exp_list:
                rec[inf] = {
                    'geomean': '',
                    'num_ref_name': num_ref_name,
                }
                csv_results.append({
                    'infection': inf,
                    'test': test,
                    'geomean': '',
                    'num_ref_name': num_ref_name
                })
                continue
            geomean = gmean(
                titer_list,
                weights=num_exp_list
            )
            rec[inf] = {
                    'geomean': round_number(geomean),
                    'num_ref_name': num_ref_name,
                }
            csv_results.append({
                'infection': inf,
                'test': test,
                'geomean': round_number(geomean),
                'num_ref_name': num_ref_name
            })

        json_results.append(rec)

    dump_csv(dst, csv_results)

    dst = dst.with_suffix('.json').with_stem('table_cross_neutral')
    dst = dst.parent.parent / dst.name
    dump_json(dst, json_results)
