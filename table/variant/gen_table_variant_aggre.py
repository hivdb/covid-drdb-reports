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


SQL = """
SELECT
    s.ref_name,
    s.variant_name,
    s.fold
FROM
    susc_results AS s,
    {rxtype} as rx
ON
    s.ref_name = rx.ref_name
    AND s.rx_name = rx.rx_name
WHERE
    s.cumulative_count > 1
"""

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
        sql = SQL.format(rxtype=rxtype)

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

        result.append(get_fold_results(single_mut_group, 'single', plasma))
        result.append(get_fold_results(combo_mut_group, 'combo', plasma))

    result.sort(key=itemgetter(
        'plasma',
        'mut_type',
        ))
    dump_csv(save_path, result)


def get_fold_results(mut_group, mut_type, plasma):
    fold_result = []
    for variant, r_list in mut_group.items():
        fold_result += [
            (r['ref_name'], r['fold']) for r in r_list if r['fold']]

    num_result = len(fold_result)
    num_ref = len(set(r[0] for r in fold_result))
    fold = [r[1] for r in fold_result]
    median_fold = median(fold) if fold else ''

    return {
        'mut_type': mut_type,
        'plasma': plasma,
        'num_ref': num_ref,
        'num_result': num_result,
        'median_fold': median_fold
    }
