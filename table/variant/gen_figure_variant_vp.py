from collections import defaultdict
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from plasma.preset import RX_VP
from .preset import ONE_MUT_VARIANT
from .preset import COMBO_MUT_VARIANT
from variant.preset import CONTROL_VARIANTS_SQL

SQL = """
SELECT
    r.vaccine_name,
    r.vaccine_type,
    s.iso_name,
    s.fold,
    s.cumulative_count as count
FROM
    susc_results as s
INNER JOIN {rx_vaccine} as r ON
    s.ref_name = r.ref_name
    AND s.rx_name = r.rx_name
WHERE
    s.potency_type IN ('IC50', 'NT50')
    AND
    s.control_iso_name in {control_variants}
    AND
    r.vaccine_name IS NOT NULL
    AND
    s.fold IS NOT NULL;
""".format(
    control_variants=CONTROL_VARIANTS_SQL,
    rx_vaccine=RX_VP,
)


def gen_figure_variant_vp(conn):
    by_variant(
        conn,
        'indiv',
        save_path=DATA_FILE_PATH / 'figure_variant_indiv_vp.csv'
    )
    by_variant(
        conn,
        'combo',
        save_path=DATA_FILE_PATH / 'figure_variant_combo_vp.csv'
    )


def by_variant(conn, indiv_or_combo, save_path):
    if indiv_or_combo == 'indiv':
        variant_mapper = ONE_MUT_VARIANT
    else:
        variant_mapper = COMBO_MUT_VARIANT
    cursor = conn.cursor()

    cursor.execute(SQL)

    variant_group = defaultdict(list)
    for rec in cursor.fetchall():
        variant = rec['iso_name']
        variant = variant_mapper.get(variant)
        if not variant:
            continue
        variant = variant['disp']
        variant_group[variant].append(rec)

    record_list = []
    for variant, rlist in variant_group.items():
        vacc_group = defaultdict(list)
        for rec in rlist:
            vacc_name = rec['vaccine_name']
            vacc_group[vacc_name].append(rec)

        for vacc_name, rx_list in vacc_group.items():

            if indiv_or_combo == 'indiv':
                variant_info = ONE_MUT_VARIANT.get(variant)
                for r in rx_list:
                    record_list.append({
                        'variant': variant,
                        'vaccine': vacc_name,
                        'vaccine_type': r['vaccine_type'],
                        'RefAA': variant_info['ref_aa'],
                        'Position': variant_info['position'],
                        'AA': variant_info['aa'],
                        'Domain': variant_info['domain'],
                        'fold': r['fold'],
                    })
            else:
                variant_info = COMBO_MUT_VARIANT.get(variant)
                nickname = variant_info['nickname']
                for r in rx_list:
                    record_list.append({
                        'variant': variant,
                        'nickname': nickname,
                        'vaccine': vacc_name,
                        'vaccine_type': r['vaccine_type'],
                        'fold': r['fold'],
                    })

    record_list.sort(key=itemgetter(
        'variant',
        'vaccine',
        ))

    dump_csv(save_path, record_list)
