from collections import defaultdict
from statistics import median
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import INDIV_VARIANT
from .preset import MULTI_VARIANT


SQL = """
SELECT
    r.vaccine_name,
    s.variant_name,
    s.fold,
    s.cumulative_count as count
FROM
    susc_results as s
INNER JOIN rx_immu_plasma as r ON
    s.ref_name = r.ref_name
    AND s.rx_name = r.rx_name
"""


def gen_table_variant_vp(conn):
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
        vacc_group = defaultdict(list)
        for rec in rlist:
            vacc_name = rec['vaccine_name']
            vacc_group[vacc_name].append(rec)

        for vacc_name, rx_list in vacc_group.items():
            all_fold = [
                [r['fold']] * r['count']
                for r in rx_list
                ]
            all_fold = [i for j in all_fold for i in j if i]
            median_fold = median(all_fold) if all_fold else ''

            num_results = sum([r['count'] for r in rx_list] + [0])

            record_list.append({
                'variant': variant,
                'vaccine': vacc_name,
                'median_fold': median_fold,
                'samples': num_results,
            })

    record_list.sort(key=itemgetter(
        'variant',
        'vaccine',
        ))

    save_path = DATA_FILE_PATH / 'table_variant_vp_figure.csv'
    dump_csv(save_path, record_list)
