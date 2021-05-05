from collections import defaultdict
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import INDIV_VARIANT
from .preset import COMBO_VARIANT
from variant.preset import CONTROL_VARIANTS_SQL

SQL = """
SELECT
    s.variant_name,
    s.fold,
    s.cumulative_count as count
FROM
    susc_results as s
INNER JOIN rx_conv_plasma as r ON
    s.ref_name = r.ref_name
    AND s.rx_name = r.rx_name
WHERE
    s.inhibition_pcnt != 90
    AND
    s.control_variant_name in {control_variants}
    AND
    s.fold IS NOT NULL;
""".format(control_variants=CONTROL_VARIANTS_SQL)


def gen_figure_variant_cp(conn):
    by_variant(
        conn,
        'indiv',
        save_path=DATA_FILE_PATH / 'figure_variant_indiv_cp.csv'
    )
    by_variant(
        conn,
        'combo',
        save_path=DATA_FILE_PATH / 'figure_variant_combo_cp.csv'
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
    for variant, rx_list in variant_group.items():

        if indiv_or_combo == 'indiv':
            variant_info = INDIV_VARIANT.get(variant)
            for r in rx_list:
                record_list.append({
                    'variant': variant,
                    'RefAA': variant_info['ref_aa'],
                    'Position': variant_info['position'],
                    'AA': variant_info['aa'],
                    'Domain': variant_info['domain'],
                    'fold': r['fold'],
                })
        else:
            variant_info = COMBO_VARIANT.get(variant)
            nickname = variant_info['nickname']
            for r in rx_list:
                record_list.append({
                    'variant': variant,
                    'nickname': nickname,
                    'fold': r['fold'],
                })

    record_list.sort(key=itemgetter(
        'variant',
        ))

    dump_csv(save_path, record_list)
