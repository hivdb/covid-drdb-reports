from collections import defaultdict
from statistics import median
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from mab.preset import MAB_RENAME
from .preset import INDIV_VARIANT
from .preset import MULTI_VARIANT
from mab.preset import ANTIBODY_TARGET_SQL
from preset import round_number


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
    by_variant(
        conn,
        'indiv',
        save_path=DATA_FILE_PATH / 'table_variant_indiv_mab_figure.csv'
    )
    by_variant(
        conn,
        'multi',
        save_path=DATA_FILE_PATH / 'table_variant_multi_mab_figure.csv'
    )


def by_variant(conn, indiv_or_multi, save_path):
    if indiv_or_multi == 'indiv':
        variant_mapper = INDIV_VARIANT
    else:
        variant_mapper = MULTI_VARIANT

    cursor = conn.cursor()
    cursor.execute(SQL)

    variant_group = defaultdict(list)
    for rec in cursor.fetchall():
        variant = rec['variant_name']
        variant = variant_mapper.get(variant, )
        if not variant:
            continue
        variant = variant['disp']
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
            target = rx_list[0]['target']

            all_fold = [
                [r['fold']] * r['count']
                for r in rx_list
                ]
            all_fold = [i for j in all_fold for i in j if i]
            median_fold = round_number(median(all_fold)) if all_fold else ''

            num_results = sum([r['count'] for r in rx_list] + [0])

            if indiv_or_multi == 'indiv':
                variant_info = INDIV_VARIANT.get(variant)
                record_list.append({
                    'variant': variant,
                    'RefAA': variant_info['ref_aa'],
                    'Position': variant_info['position'],
                    'AA': variant_info['aa'],
                    'Domain': variant_info['domain'],
                    'mab': rx_name,
                    'avail': avail,
                    'target': target,
                    'median_fold': median_fold,
                    'samples': num_results,
                })
            else:
                record_list.append({
                    'variant': variant,
                    'mab': rx_name,
                    'avail': avail,
                    'target': target,
                    'median_fold': median_fold,
                    'samples': num_results,
                })

    record_list.sort(key=itemgetter(
        'variant',
        'mab',
        'avail',
        ))
    dump_csv(save_path, record_list)
