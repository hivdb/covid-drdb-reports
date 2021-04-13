from collections import defaultdict
from statistics import median
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import INDIV_VARIANT
from .preset import MULTI_VARIANT

from mab.preset import MAB_RENAME
from mab.preset import RX_MAB

from preset import round_number

from resistancy import is_partial_resistant
from resistancy import is_resistant
from resistancy import is_susc

SQL = """
SELECT
    rx.ab_name,
    s.variant_name,
    s.fold,
    s.cumulative_count as count,
    rx.availability as avail,
    rx.pdb_id as pdb,
    rx.target as target
FROM
    susc_results as s,
    ({rx_type}) as rx
ON
    s.ref_name = rx.ref_name
    AND s.rx_name = rx.rx_name
WHERE
    s.control_variant_name IN ('Control', 'Wuhan', 'S:614G');
""".format(rx_type=RX_MAB)


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
        variant = variant_mapper.get(variant)
        if not variant:
            continue
        variant = variant['disp']
        variant_group[variant].append(rec)

    record_list = []
    for variant, rlist in variant_group.items():
        rx_group = defaultdict(list)
        for rec in rlist:
            ab_name = rec['ab_name']
            rx_group[ab_name].append(rec)

        for ab_name, rx_list in rx_group.items():
            ab_name = MAB_RENAME.get(ab_name, ab_name)
            avail = rx_list[0]['avail']
            target = rx_list[0]['target']

            all_fold = [(r['fold'], r['count']) for r in rx_list if r['fold']]
            num_s = sum([r[1] for r in all_fold if is_susc(r[0])])
            num_i = sum([r[1] for r in all_fold if is_partial_resistant(r[0])])
            num_r = sum([r[1] for r in all_fold if is_resistant(r[0])])

            all_fold = [[i[0]] * i[1] for i in all_fold]
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
                    'mab': ab_name,
                    'avail': avail,
                    'target': target,
                    'median_fold': median_fold,
                    'samples': num_results,
                    'S': num_s,
                    'I': num_i,
                    'R': num_r,
                })
            else:
                record_list.append({
                    'variant': variant,
                    'mab': ab_name,
                    'avail': avail,
                    'target': target,
                    'median_fold': median_fold,
                    'samples': num_results,
                    'S': num_s,
                    'I': num_i,
                    'R': num_r,
                })

    record_list.sort(key=itemgetter(
        'variant',
        'mab',
        'avail',
        ))
    dump_csv(save_path, record_list)
