from preset import DATA_FILE_PATH
from preset import load_csv
from preset import dump_csv
from .preset import IGNORE_VACCINE_NAME

from .common import get_sample_number_pair
from collections import defaultdict
from itertools import product
from variant.preset import SINGLE_S_MUTATION_ISOLATES


CP_TITER_VARIANT_SQL = """
SELECT
    c.var_name iso_name,
    "CP" rx_name,
    a.ref_name,
    a.potency titer,
    b.timing month,
    a.cumulative_count num_result
FROM
    rx_potency a,
    rx_conv_plasma b,
    isolates c
WHERE
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    a.iso_name = c.iso_name
    AND
    c.var_name IS NOT NULL
    AND
    a.potency_type = 'NT50'
"""

CP_TITER_SINGLE_MUT_SQL = """
SELECT
    c.single_mut_name iso_name,
    "CP" rx_name,
    a.ref_name,
    a.potency titer,
    b.timing month,
    a.cumulative_count num_result
FROM
    rx_potency a,
    rx_conv_plasma b,
    ({single_s_mutation_isolates}) c
WHERE
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    a.iso_name = c.iso_name
    AND
    a.potency_type = 'NT50'
"""

VP_TITER_VARIANT_SQL = """
SELECT
    c.var_name iso_name,
    vaccine_name rx_name,
    a.ref_name,
    a.potency titer,
    b.timing month,
    a.cumulative_count num_result
FROM
    rx_potency a,
    rx_vacc_plasma b,
    isolates c
WHERE
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    a.iso_name = c.iso_name
    AND
    c.var_name IS NOT NULL
    AND
    a.potency_type = 'NT50'
"""

VP_TITER_SINGLE_MUT_SQL = """
SELECT
    c.single_mut_name iso_name,
    vaccine_name as rx_name,
    a.ref_name,
    a.potency titer,
    b.timing month,
    a.cumulative_count num_result
FROM
    rx_potency a,
    rx_vacc_plasma b,
    ({single_s_mutation_isolates}) c
WHERE
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    a.iso_name = c.iso_name
    AND
    a.potency_type = 'NT50'
"""


ISO_NAME_LIST = [
    "Alpha",
    "Beta",
    "Gamma",
    "Delta",
    "N501Y",
    "E484K",
    "L452R"
]

RX_NAME_LIST = [
    "CP",
    "BNT162b2",
    'mRNA-1273',
    'AZD1222',
    'Ad26.COV2.S',
    'NVX-CoV2373',
    'BBV152',
    'CoronaVac',
    'BBIBP-CorV',
    'Sputnik V',
    'MVC-COV1901',
    'ZF2001',
]

LEVEL = [
    'low',
    'middle',
    'high',
]

MONTH = [
    '1',
    '2-6',
    '>6',
]


def gen_figure_plasma_titer(
        conn,
        save_path=DATA_FILE_PATH / 'figure_plasma_titer.csv'):
    sql_tmpl = " UNION ALL ".join([
        CP_TITER_VARIANT_SQL,
        CP_TITER_SINGLE_MUT_SQL,
        VP_TITER_VARIANT_SQL,
        VP_TITER_SINGLE_MUT_SQL
    ])

    records = []
    used_header_combo = []

    sql = sql_tmpl.format(
        single_s_mutation_isolates=SINGLE_S_MUTATION_ISOLATES
        )

    cursor = conn.cursor()
    cursor.execute(sql)

    rows = []
    for row in cursor.fetchall():
        rec = {}
        for key in row.keys():
            rec[key] = row[key]
        rows.append(rec)

    iso_name_group = defaultdict(list)
    for rec in rows:
        iso_name = rec['iso_name']
        iso_name_group[iso_name].append(rec)

    for iso_name, iso_rec_list in iso_name_group.items():

        rx_group = defaultdict(list)
        for rec in iso_rec_list:
            rx_name = rec['rx_name']
            rx_group[rx_name].append(rec)

        for rx_name, rx_rec_list in rx_group.items():

            month_group = defaultdict(list)
            for rec in rx_rec_list:
                month = rec['month']
                if month < 2:
                    month = '1'
                elif month >= 2 and month <= 6:
                    month = '2-6'
                else:
                    month = '>6'
                month_group[month].append(rec)

            for month, month_rec_list in month_group.items():
                low = [r for r in month_rec_list if int(r['titer']) <= 40]
                middle = [
                    r for r in month_rec_list if int(r['titer']) > 40
                    and int(r['titer']) <= 100]
                high = [r for r in month_rec_list if int(r['titer']) > 100]

                records.append({
                    'iso_name': iso_name,
                    'rx_name': rx_name,
                    'level': 'low',
                    'month': month,
                    'num_study': len(
                        set(r['ref_name'] for r in month_rec_list)),
                    'num_sample': sum([r['num_result'] for r in low])
                })

                used_header_combo.append((
                    iso_name,
                    rx_name,
                    'low',
                    month
                ))

                records.append({
                    'iso_name': iso_name,
                    'rx_name': rx_name,
                    'level': 'middle',
                    'month': month,
                    'num_study': len(
                        set(r['ref_name'] for r in month_rec_list)),
                    'num_sample': sum([r['num_result'] for r in middle])
                })

                used_header_combo.append((
                    iso_name,
                    rx_name,
                    'middle',
                    month
                ))

                records.append({
                    'iso_name': iso_name,
                    'rx_name': rx_name,
                    'level': 'high',
                    'month': month,
                    'num_study': len(
                        set(r['ref_name'] for r in month_rec_list)),
                    'num_sample': sum([r['num_result'] for r in high])
                })

                used_header_combo.append((
                    iso_name,
                    rx_name,
                    'high',
                    month
                ))

    for store_key in product(
            ISO_NAME_LIST, RX_NAME_LIST, LEVEL, MONTH):
        if store_key in used_header_combo:
            continue
        used_header_combo.append(store_key)
        iso_name, rx_name, level, month = store_key
        records.append({
                'iso_name': iso_name,
                'rx_name': rx_name,
                'level': level,
                'month': month,
                'num_study': 0,
                'num_sample': 0
            })

    dump_csv(save_path, records)
