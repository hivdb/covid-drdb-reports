from preset import DATA_FILE_PATH
from preset import load_csv
from preset import dump_csv
from .preset import IGNORE_VACCINE_NAME

from .common import get_sample_number_pair
from collections import defaultdict
from variant.preset import SINGLE_S_MUTATION_ISOLATES


CP_TITER_VARIANT_SQL = """
SELECT
    c.var_name iso_name,
    "CP" rx_name,
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
""".format(
    single_s_mutation_isolates=SINGLE_S_MUTATION_ISOLATES
)

VP_TITER_VARIANT_SQL = """
SELECT
    c.var_name iso_name,
    vaccine_name rx_name,
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
""".format(
    single_s_mutation_isolates=SINGLE_S_MUTATION_ISOLATES
)


def gen_figure_plasma_titer(
        conn,
        save_path=DATA_FILE_PATH / 'figure_plasma_titer.csv'):
    sql = " UNION ".join([
        CP_TITER_VARIANT_SQL,
        CP_TITER_SINGLE_MUT_SQL,
        VP_TITER_VARIANT_SQL,
        VP_TITER_SINGLE_MUT_SQL
    ])

    cursor = conn.cursor()
    cursor.execute(sql)

    records = []
    for row in cursor.fetchall():
        rec = {}
        for key in row.keys():
            rec[key] = row[key]
        records.append(rec)

    results = []

    iso_group = defaultdict(list)
    for rec in records:
        iso_name = rec['iso_name']
        iso_group[iso_name].append(rec)

    for iso, rec_list in iso_group.items():

        rx_group = defaultdict(list)
        for rec in rec_list:
            rx_name = rec['rx_name']
            rx_group[rx_name].append(rec)

        for rx_name, rx_rec_list in rx_group.items():
            low = [r for r in rx_rec_list if int(r['titer']) <= 40]
            middle = [r for r in rx_rec_list if int(r['titer']) > 40 and int(r['titer']) <= 100]
            high = [r for r in rx_rec_list if int(r['titer']) > 100]

            results.append({
                'iso_name': iso,
                'rx_name': rx_name,
                'level': 'l',
                'count': sum([r['num_result'] for r in low])
            })

            results.append({
                'iso_name': iso,
                'rx_name': rx_name,
                'level': 'm',
                'count': sum([r['num_result'] for r in middle])
            })

            results.append({
                'iso_name': iso,
                'rx_name': rx_name,
                'level': 'h',
                'count': sum([r['num_result'] for r in high])
            })

    dump_csv(save_path, results)
