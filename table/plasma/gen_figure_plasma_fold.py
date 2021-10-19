from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import group_records_by
from resistancy import get_susceptibility
from collections import defaultdict
from itertools import product
from .preset import INFECTED_VACCINEE


CP_FOLD_VARIANT_SQL = """
SELECT
    c.var_name iso_name,
    "CP" rx_name,
    a.ref_name,
    a.fold_cmp || a.fold fold,
    b.timing month,
    'infected' infection,
    a.cumulative_count num_result
FROM
    susc_results_view a,
    rx_conv_plasma b,
    isolates c,
    subjects d
WHERE
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    a.iso_name = c.iso_name
    AND
    c.var_name IS NOT NULL
    AND
    b.ref_name = d.ref_name
    AND
    b.subject_name = d.subject_name
    AND
    d.subject_species = 'Human'
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

VP_FOLD_VARIANT_SQL = """
SELECT
    c.var_name iso_name,
    b.vaccine_name rx_name,
    a.ref_name,
    a.fold_cmp || a.fold fold,
    b.timing month,
    'naive' infection,
    a.cumulative_count num_result
FROM
    susc_results_view a,
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

VP_FOLD_INFECTED_VARIANT_SQL = """
SELECT
    c.var_name iso_name,
    b.vaccine_name rx_name,
    a.ref_name,
    a.fold_cmp || a.fold fold,
    b.timing month,
    'infected' infection,
    a.cumulative_count num_result
FROM
    susc_results_view a,
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


def gen_figure_plasma_fold(
        conn,
        save_path=DATA_FILE_PATH / 'figure' / 'figure_plasma_fold.csv',
        save_path_indiv=DATA_FILE_PATH / 'figure' / 'figure_plasma_fold_indiv.csv'):
    sql_tmpl = " UNION ALL ".join([
        CP_FOLD_VARIANT_SQL,
        VP_FOLD_VARIANT_SQL,
        VP_FOLD_INFECTED_VARIANT_SQL,
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

    padding_records(records, used_header_combo)

    dump_csv(save_path, records)

    records, used_header_combo = _add_records(rows, True)

    padding_records(records, used_header_combo)

    dump_csv(save_path_indiv, records)


def _add_records(records, indiv_data_only=False):
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

        susc = [
            r for r in rec_list
            if get_susceptibility(r['fold']) == 'susceptible']
        partial = [
            r for r in rec_list
            if get_susceptibility(r['fold']) == 'partial-resistance']
        resist = [
            r for r in rec_list
            if get_susceptibility(r['fold']) == 'resistant']

        if indiv_data_only:
            num_result = sum([
                r['num_result'] for r in susc if r['num_result'] == 1])
            num_study = len(
                set(
                    r['ref_name'] for r in rec_list if r['num_result'] == 1))
        else:
            num_result = sum([r['num_result'] for r in susc])
            num_study = len(
                set(r['ref_name'] for r in rec_list))

        data_points.append({
            'variant': iso_name,
            'rx_name': rx_name,
            'susc': 'S',
            'month': month,
            'infection': infection,
            'num_study': num_study,
            'num_result': num_result
        })

        used_header_combo.append((
            iso_name,
            rx_name,
            'S',
            month,
            infection
        ))

        if indiv_data_only:
            num_result = sum([
                r['num_result'] for r in partial if r['num_result'] == 1])
            num_study = len(
                set(
                    r['ref_name'] for r in rec_list if r['num_result'] == 1))
        else:
            num_result = sum([r['num_result'] for r in partial])
            num_study = len(
                set(r['ref_name'] for r in rec_list))

        data_points.append({
            'variant': iso_name,
            'rx_name': rx_name,
            'susc': 'I',
            'month': month,
            'infection': infection,
            'num_study': num_study,
            'num_result': num_result
        })

        used_header_combo.append((
            iso_name,
            rx_name,
            'I',
            month,
            infection
        ))

        if indiv_data_only:
            num_result = sum([
                r['num_result'] for r in resist if r['num_result'] == 1])
            num_study = len(
                set(
                    r['ref_name'] for r in rec_list if r['num_result'] == 1))
        else:
            num_result = sum([r['num_result'] for r in resist])
            num_study = len(
                set(r['ref_name'] for r in rec_list))

        data_points.append({
            'variant': iso_name,
            'rx_name': rx_name,
            'susc': 'R',
            'month': month,
            'infection': infection,
            'num_study': num_study,
            'num_result': num_result
        })

        used_header_combo.append((
            iso_name,
            rx_name,
            'R',
            month,
            infection
        ))

    return data_points, used_header_combo


def padding_records(records, used_header_combo):
    for store_key in product(
                ISO_NAME_LIST,
                RX_NAME_LIST,
                FOLD_LEVEL,
                MONTH,
                INFECTION):
        if store_key in used_header_combo:
            continue
        used_header_combo.append(store_key)
        iso_name, rx_name, level, month, infection = store_key
        records.append({
                'variant': iso_name,
                'rx_name': rx_name,
                'susc': level,
                'month': month,
                'infection': infection,
                'num_study': 0,
                'num_result': 0
            })
