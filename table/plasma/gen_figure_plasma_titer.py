from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import group_records_by
from collections import defaultdict
from itertools import product
from .preset import INFECTED_VACCINEE


CP_TITER_VARIANT_SQL = """
SELECT
    c.var_name iso_name,
    "CP" rx_name,
    a.ref_name,
    a.potency titer,
    b.timing month,
    'infected' infection,
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
    AND
    NOT EXISTS (
        SELECT
            1
        FROM ({infected_vaccinee}) d
        WHERE
            b.ref_name = d.ref_name
            AND
            b.subject_name = d.subject_name
    )
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
    b.vaccine_name rx_name,
    a.ref_name,
    a.potency titer,
    b.timing month,
    'naive' infection,
    a.cumulative_count num_result
FROM
    rx_potency a,
    rx_vacc_plasma b,
    isolates c,
    vaccines d
WHERE
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    a.iso_name = c.iso_name
    AND
    c.var_name IS NOT NULL
    AND
    b.vaccine_name = d.vaccine_name
    AND
    b.dosage = d.st_shot
    AND
    a.potency_type = 'NT50'
    AND
    NOT EXISTS (
        SELECT
            1
        FROM ({infected_vaccinee}) f
        WHERE
            b.ref_name = f.ref_name
            AND
            b.subject_name = f.subject_name
    )
"""

VP_TITER_INFECTED_VARIANT_SQL = """
SELECT
    c.var_name iso_name,
    b.vaccine_name rx_name,
    a.ref_name,
    a.potency titer,
    b.timing month,
    'infected' infection,
    a.cumulative_count num_result
FROM
    rx_potency a,
    rx_vacc_plasma b,
    isolates c,
    vaccines d
WHERE
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    a.iso_name = c.iso_name
    AND
    c.var_name IS NOT NULL
    AND
    b.vaccine_name = d.vaccine_name
    AND
    b.dosage = d.st_shot
    AND
    a.potency_type = 'NT50'
    AND
    EXISTS (
        SELECT
            1
        FROM ({infected_vaccinee}) f
        WHERE
            b.ref_name = f.ref_name
            AND
            b.subject_name = f.subject_name
    )
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
    '>=2',
]

INFECTION = [
    'infected',
    'naive'
]


def gen_figure_plasma_titer(
        conn,
        save_path=DATA_FILE_PATH / 'figure_plasma_titer.csv'):
    sql_tmpl = " UNION ALL ".join([
        CP_TITER_VARIANT_SQL,
        VP_TITER_VARIANT_SQL,
        VP_TITER_INFECTED_VARIANT_SQL,
    ])

    sql = sql_tmpl.format(infected_vaccinee=INFECTED_VACCINEE)

    cursor = conn.cursor()
    cursor.execute(sql)

    rows = []
    for row in cursor.fetchall():
        rec = {}
        for key in row.keys():
            rec[key] = row[key]
        rows.append(rec)

    records, used_header_combo = _add_records(rows)

    for store_key in product(
            ISO_NAME_LIST,
            RX_NAME_LIST,
            LEVEL,
            MONTH,
            INFECTION):
        if store_key in used_header_combo:
            continue
        used_header_combo.append(store_key)
        iso_name, rx_name, level, month, infection = store_key
        records.append({
                'variant': iso_name,
                'rx_name': rx_name,
                'level': level,
                'month': month,
                'infection': infection,
                'num_study': 0,
                'num_result': 0
            })

    dump_csv(save_path, records)

    indiv_records = []
    used_key = []
    for i in rows:
        cumu = i['num_result']
        i['num_result'] = 1
        del i['ref_name']
        month = i['month']
        if month < 2:
            month = '1'
        else:
            month = '>=2'
        i['month'] = month

        for _ in range(int(cumu)):
            indiv_records.append(i)

        key = i['iso_name'], i['rx_name'], i['month'], i['infection']
        if key in used_key:
            continue
        used_key.append(key)

    for key in product(ISO_NAME_LIST, RX_NAME_LIST, MONTH, INFECTION):
        if key in used_key:
            continue
        used_key.append(key)
        iso_name, rx_name, month, infection = key

        indiv_records.append({
            'iso_name': iso_name,
            'rx_name': rx_name,
            'titer': '',
            'month': month,
            'infection': infection,
            'num_result': 1
        })

    dump_csv(DATA_FILE_PATH / 'figure_plasma_titer_points.csv', indiv_records)


def _add_records(records):
    data_points = []
    used_header_combo = []

    grouped_records = group_records_by(records, 'iso_name')
    new_grouped_records = {}

    for key in ['rx_name', 'infection']:
        for group_key, rec_list in grouped_records.items():
            if type(group_key) == str:
                group_key = [group_key]

            sub_group = group_records_by(rec_list, key)
            for sub_key, sub_rec_list in sub_group.items():
                new_group_key = tuple(list(group_key) + [sub_key])
                new_grouped_records[new_group_key] = sub_rec_list

        grouped_records = new_grouped_records
        new_grouped_records = {}

    for group_key, rec_list in grouped_records.items():
        month_group = defaultdict(list)
        for rec in rec_list:
            month = rec['month']
            if month < 2:
                month = '1'
            else:
                month = '>=2'
            month_group[month].append(rec)

        for sub_key, sub_rec_list in month_group.items():
            new_group_key = tuple(list(group_key) + [sub_key])
            new_grouped_records[new_group_key] = sub_rec_list
    grouped_records = new_grouped_records
    new_grouped_records = {}

    for group_key, rec_list in grouped_records.items():
        iso_name, rx_name, infection, month = group_key

        low = [
            r for r in rec_list
            if int(r['titer']) <= 40]
        middle = [
            r for r in rec_list
            if int(r['titer']) > 40
            and int(r['titer']) <= 100]
        high = [
            r for r in rec_list
            if int(r['titer']) > 100]

        data_points.append({
            'variant': iso_name,
            'rx_name': rx_name,
            'level': 'low',
            'month': month,
            'infection': infection,
            'num_study': len(
                set(r['ref_name'] for r in rec_list)),
            'num_result': sum([r['num_result'] for r in low])
        })

        used_header_combo.append((
            iso_name,
            rx_name,
            'low',
            month,
            infection
        ))

        data_points.append({
            'variant': iso_name,
            'rx_name': rx_name,
            'level': 'middle',
            'month': month,
            'infection': infection,
            'num_study': len(
                set(r['ref_name'] for r in rec_list)),
            'num_result': sum([r['num_result'] for r in middle])
        })

        used_header_combo.append((
            iso_name,
            rx_name,
            'middle',
            month,
            infection
        ))

        data_points.append({
            'variant': iso_name,
            'rx_name': rx_name,
            'level': 'high',
            'month': month,
            'infection': infection,
            'num_study': len(
                set(r['ref_name'] for r in rec_list)),
            'num_result': sum([r['num_result'] for r in high])
        })

        used_header_combo.append((
            iso_name,
            rx_name,
            'high',
            month,
            infection
        ))

    return data_points, used_header_combo
