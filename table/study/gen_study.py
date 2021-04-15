from preset import DATA_FILE_PATH
from preset import dump_csv
from collections import defaultdict

from variant.preset import INDIV_VARIANT
from variant.preset import COMBO_VARIANT
from variant.preset import group_by_variant
from mab.preset import RX_MAB
from variant.preset import CONTROL_VARIANTS_SQL

SQL = """
SELECT
    rx.ref_name,
    s.variant_name
FROM
    susc_results as s,
    {rx_type} as rx
ON
    s.ref_name = rx.ref_name
    AND s.rx_name = rx.rx_name
WHERE
    s.control_variant_name IN {control_variants}
"""

# MAB_SQL = """
# SELECT
#     rx.ref_name,
#     s.variant_name,
#     rx.target
# FROM
#     susc_results as s,
#     {rx_type} as rx,
# ON
#     s.ref_name = rx.ref_name
#     AND s.rx_name = rx.rx_name
# WHERE
#     s.control_variant_name IN {control_variants}
# """


RX_TYPE = {
    'rx_conv_plasma': DATA_FILE_PATH / 'summary_study_cp.csv',
    'rx_immu_plasma': DATA_FILE_PATH / 'summary_study_vp.csv',
    '({})'.format(RX_MAB): DATA_FILE_PATH / 'summary_study_mab.csv',
}


def gen_study(conn):

    cursor = conn.cursor()

    for rx_type, save_path in RX_TYPE.items():
        sql = SQL.format(
            rx_type=rx_type,
            control_variants=CONTROL_VARIANTS_SQL,
            )

        cursor.execute(sql)

        records = cursor.fetchall()

        indiv_records, combo_records = group_by_variant(records)

        results = []
        results.extend(get_indiv_studies(indiv_records))
        results.extend(get_combo_studies(combo_records))

        dump_csv(save_path, results)


def get_indiv_studies(records):
    results = []

    domain_groups = defaultdict(list)
    for variant, rec in records.items():
        variant_info = INDIV_VARIANT.get(variant)
        domain = variant_info['domain']
        domain_groups[domain].extend(rec)

    for domain, rx_list in domain_groups.items():
        studies = [r['ref_name'] for r in rx_list]
        num_studies = len(set(studies))

        results.append({
            'variant_type': 'indiv',
            'domain': domain,
            'studies': num_studies,
        })

    return results


def get_combo_studies(records):
    results = []

    studies = []

    for _, rx_list in records.items():
        studies.extend([r['ref_name'] for r in rx_list])

    results.append({
        'variant_type': 'combo',
        'domain': '',
        'studies': len(set(studies))
    })

    return results
