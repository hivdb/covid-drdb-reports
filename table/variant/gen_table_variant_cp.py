from collections import defaultdict
from statistics import median
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import INDIV_VARIANT
from .preset import MULTI_VARIANT
from preset import round_number

SQL = """
SELECT
    s.variant_name,
    s.fold,
    s.cumulative_count as count
FROM
    susc_results as s
INNER JOIN rx_immu_plasma as r ON
    s.ref_name = r.ref_name
    AND s.rx_name = r.rx_name
WHERE s.control_variant_name in ('Wuhan', 'Control')
"""


def gen_table_variant_cp(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)

    variant_group = defaultdict(list)
    for rec in cursor.fetchall():
        variant = rec['variant_name']
        variant = INDIV_VARIANT.get(variant, MULTI_VARIANT.get(variant))
        if not variant:
            continue
        variant_group[variant].append(rec)

    record_list = []
    for variant, rlist in variant_group.items():
        all_fold = [
            [r['fold']] * r['count']
            for r in rlist
            ]
        all_fold = [i for j in all_fold for i in j if i]
        median_fold = round_number(median(all_fold)) if all_fold else ''

        num_results = sum([r['count'] for r in rlist] + [0])

        record_list.append({
            'variant': variant,
            'median_fold': median_fold,
            'samples': num_results,
        })

    record_list.sort(key=itemgetter(
        'variant',
        ))

    save_path = DATA_FILE_PATH / 'table_variant_cp_figure.csv'
    dump_csv(save_path, record_list)
