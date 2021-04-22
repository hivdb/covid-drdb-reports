from collections import defaultdict
from statistics import median
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import INDIV_VARIANT
from .preset import COMBO_VARIANT

from mab.preset import MAB_RENAME
from mab.preset import RX_MAB

from resistancy import round_fold

from resistancy import is_partial_resistant
from resistancy import is_resistant
from resistancy import is_susc
from variant.preset import CONTROL_VARIANTS_SQL


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
    s.control_variant_name IN {control_variants}
    AND s.fold IS NOT NULL;
""".format(rx_type=RX_MAB, control_variants=CONTROL_VARIANTS_SQL)


def gen_figure_variant_mab(conn):
    by_variant(
        conn,
        'indiv',
        save_path=DATA_FILE_PATH / 'figure_variant_indiv_mab.csv'
    )
    by_variant(
        conn,
        'combo',
        save_path=DATA_FILE_PATH / 'figure_variant_combo_mab.csv'
    )


def by_variant(conn, indiv_or_combo, save_path):
    if indiv_or_combo == 'indiv':
        variant_mapper = INDIV_VARIANT
    else:
        variant_mapper = COMBO_VARIANT

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

            if indiv_or_combo == 'indiv':
                variant_info = INDIV_VARIANT.get(variant)
                for r in rx_list:
                    record_list.append({
                        'variant': variant,
                        'RefAA': variant_info['ref_aa'],
                        'Position': variant_info['position'],
                        'AA': variant_info['aa'],
                        'Domain': variant_info['domain'],
                        'mab': ab_name,
                        'avail': avail,
                        'target': target,
                        'fold': r['fold'],
                    })
            else:
                variant_info = COMBO_VARIANT.get(variant)
                nickname = variant_info['nickname']
                for r in rx_list:
                    record_list.append({
                        'variant': variant,
                        'nickname': nickname,
                        'mab': ab_name,
                        'avail': avail,
                        'target': target,
                        'fold': r['fold'],
                    })

    record_list.sort(key=itemgetter(
        'variant',
        'mab',
        'avail',
        ))
    dump_csv(save_path, record_list)
