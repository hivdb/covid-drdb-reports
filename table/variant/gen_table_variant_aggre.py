from collections import defaultdict
from statistics import median
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import INDIV_VARIANT
from .preset import MULTI_VARIANT
from preset import round_number
from resistancy import is_partial_resistant
from resistancy import is_resistant
from resistancy import is_susc
from plasma.preset import AGGREGATED_RESULTS_SQL

RXTYPE = {
    'rx_immu_plasma': 'vp',
    'rx_conv_plasma': 'cp',
}


def gen_table_variant_aggre(
        conn,
        save_path=DATA_FILE_PATH / 'table_variant_aggre_plasma_figure.csv'):

    cursor = conn.cursor()

    result = []

    for rxtype, plasma in RXTYPE.items():
        sql = AGGREGATED_RESULTS_SQL.format(
            rxtype=rxtype, filters='')

        # print(sql)
        cursor.execute(sql)

        single_mut_group = defaultdict(list)
        combo_mut_group = defaultdict(list)
        for rec in cursor.fetchall():
            variant_name = rec['variant_name']
            variant = INDIV_VARIANT.get(variant_name)
            if variant:
                variant = variant['disp']
                single_mut_group[variant].append(rec)
                continue

            variant = MULTI_VARIANT.get(variant_name)
            if not variant:
                continue
            variant = variant['disp']
            combo_mut_group[variant].append(rec)

        result += get_fold_results(single_mut_group, 'single', plasma)
        result += get_fold_results(combo_mut_group, 'combo', plasma)

    result.sort(key=itemgetter(
        'plasma',
        ))
    dump_csv(save_path, result)


def get_fold_results(mut_group, mut_type, plasma):
    results = []
    for variant, r_list in mut_group.items():
        for r in r_list:
            results.append({
                'variant': variant,
                'plasma': plasma,
                'reference': r['ref_name'],
                'fold_cmp': r['fold_cmp'],
                'median': r['fold'],
                'count': r['sample_count']
            })

    return results
