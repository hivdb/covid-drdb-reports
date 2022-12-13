from preset import DATA_FILE_PATH
from preset import dump_csv
from itertools import product
from .preset import INFECTED_VACCINEE
from .preset import RX_TITER_FOLD
from resistancy import parse_fold


CP_TITER_FOLD_VARIANT_SQL = """
SELECT
    c.var_name iso_name,
    "CP" rx_name,
    a.ref_name,
    a.fold_cmp || a.fold fold,
    a.potency titer,
    b.timing month,
    'infected' infection,
    a.cumulative_count num_result
FROM
    ({rx_titer_fold}) a,
    rx_conv_plasma b,
    isolates c,
    subjects d
WHERE
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    b.ref_name = d.ref_name
    AND
    b.subject_name = d.subject_name
    AND
    d.subject_species = 'Human'
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
        FROM ({infected_vaccinee}) inf
        WHERE
            b.ref_name = inf.ref_name
            AND
            b.subject_name = inf.subject_name
    )
"""

VP_TITER_FOLD_VARIANT_SQL = """
SELECT
    c.var_name iso_name,
    b.vaccine_name rx_name,
    a.ref_name,
    a.fold_cmp || a.fold fold,
    a.potency titer,
    b.timing month,
    'naive' infection,
    a.cumulative_count num_result
FROM
    ({rx_titer_fold}) a,
    rx_vacc_plasma b,
    isolates c,
    subjects d,
    vaccines e
WHERE
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    a.iso_name = c.iso_name
    AND
    b.ref_name = d.ref_name
    AND
    b.subject_name = d.subject_name
    AND
    d.subject_species = 'Human'
    AND
    c.var_name IS NOT NULL
    AND
    b.vaccine_name = e.vaccine_name
    AND
    b.dosage = e.st_shot
    AND
    a.potency_type = 'NT50'
    AND
    NOT EXISTS (
        SELECT
            1
        FROM ({infected_vaccinee}) inf
        WHERE
            b.ref_name = inf.ref_name
            AND
            b.subject_name = inf.subject_name
    )
"""

VP_TITER_FOLD_INFECTED_VARIANT_SQL = """
SELECT
    c.var_name iso_name,
    b.vaccine_name rx_name,
    a.ref_name,
    a.fold_cmp || a.fold fold,
    a.potency titer,
    b.timing month,
    'infected' infection,
    a.cumulative_count num_result
FROM
    ({rx_titer_fold}) a,
    rx_vacc_plasma b,
    isolates c,
    subjects d,
    vaccines e
WHERE
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    a.iso_name = c.iso_name
    AND
    b.ref_name = d.ref_name
    AND
    b.subject_name = d.subject_name
    AND
    d.subject_species = 'Human'
    AND
    c.var_name IS NOT NULL
    AND
    b.vaccine_name = e.vaccine_name
    AND
    b.dosage = e.st_shot
    AND
    a.potency_type = 'NT50'
    AND
    EXISTS (
        SELECT
            1
        FROM ({infected_vaccinee}) inf
        WHERE
            b.ref_name = inf.ref_name
            AND
            b.subject_name = inf.subject_name
    )
"""


ISO_NAME_LIST = [
    "Alpha",
    "Beta",
    "Gamma",
    "Delta",
    "Omicron/BA.1",
    "Omicron/BA.2",
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

FOLD_LEVEL = [
    'S',
    'I',
    'R',
]

MONTH = [
    '1',
    '>=2',
]

INFECTION = [
    'infected',
    'naive'
]


def gen_figure_plasma_titer_fold(
        conn,
        save_path=DATA_FILE_PATH / 'figure' / 'figure_plasma_titer_fold.csv'):
        
    raise Exception('Archived analysis')

    sql_tmpl = " UNION ALL ".join([
        CP_TITER_FOLD_VARIANT_SQL,
        VP_TITER_FOLD_VARIANT_SQL,
        VP_TITER_FOLD_INFECTED_VARIANT_SQL,
    ])


    sql = sql_tmpl.format(
        infected_vaccinee=INFECTED_VACCINEE,
        rx_titer_fold=RX_TITER_FOLD
        )

    cursor = conn.cursor()
    cursor.execute(sql)

    rows = []
    for row in cursor.fetchall():
        rec = {}
        for key in row.keys():
            rec[key] = row[key]
        rows.append(rec)

    indiv_records = []
    used_key = []
    for i in rows:
        cumu = i['num_result']

        if cumu > 1:
            continue

        i['variant'] = i['iso_name']

        month = i['month']
        if month < 2:
            month = '1'
        else:
            month = '>=2'
        i['month'] = month

        i['fold'] = parse_fold(i['fold'])
        if i['fold'] == 0:
            i['fold'] = 0.1

        for _ in range(int(cumu)):
            indiv_records.append(i)

        key = i['variant'], i['rx_name'], i['month'], i['infection']
        if key in used_key:
            continue
        used_key.append(key)

    for key in product(ISO_NAME_LIST, RX_NAME_LIST, MONTH, INFECTION):
        if key in used_key:
            continue
        used_key.append(key)
        iso_name, rx_name, month, infection = key

        indiv_records.append({
            'variant': iso_name,
            'iso_name': '',
            'ref_name': '',
            'rx_name': rx_name,
            'titer': '',
            'fold': '',
            'month': month,
            'infection': infection,
            'num_result': 1
        })

    dump_csv(
        DATA_FILE_PATH /
        'figure' / 'figure_plasma_titer_fold_points.csv', indiv_records)
