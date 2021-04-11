from collections import defaultdict
from statistics import median
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from mab.preset import MAB_RENAME
from .preset import INDIV_VARIANT
from .preset import MULTI_VARIANT
from mab.preset import ANTIBODY_TARGET_SQL


SQL = """
SELECT
    s.rx_name,
    s.variant_name,
    s.fold,
    s.cumulative_count as count,
    a.availability as avail,
    t.pdb_id as pdb,
    t.target as target
FROM
    susc_results as s
INNER JOIN rx_antibodies as r ON
    s.ref_name = r.ref_name
    AND s.rx_name = r.rx_name
LEFT JOIN antibodies as a ON a.ab_name = r.ab_name
LEFT JOIN """ + ANTIBODY_TARGET_SQL + """ as t ON a.ab_name = t.ab_name
WHERE
    r.ab_name in (
        SELECT ab_name FROM 'antibodies'
    )
"""


def gen_table_variant_mab(conn):
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
        rx_group = defaultdict(list)
        for rec in rlist:
            rx_name = rec['rx_name']
            rx_group[rx_name].append(rec)

        for rx_name, rx_list in rx_group.items():
            rx_name = MAB_RENAME.get(rx_name, rx_name)
            avail = rx_list[0]['avail']

            fold = [r['fold'] for r in rx_list]
            median_fold = median(fold)
            num_results = sum([r['count'] for r in rx_list] + [0])

            record_list.append({
                'variant': variant,
                'mab': rx_name,
                'avail': avail,
                'median': median_fold,
                'samples': num_results,
            })

    record_list.sort(key=itemgetter(
        'variant',
        'mab',
        'avail',
        ))

    save_path = DATA_FILE_PATH / 'table_variant_mab_figure.csv'
    dump_csv(save_path, record_list)
