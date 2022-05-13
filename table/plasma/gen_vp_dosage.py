from operator import itemgetter
from preset import dump_csv
from preset import dump_json
from preset import DATA_FILE_PATH
from collections import defaultdict


SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    rx.dosage,
    SUM(s.cumulative_count) num_fold
FROM
    susc_results_view s,
    rx_vacc_plasma rx
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.fold IS NOT NULL
GROUP BY
    s.ref_name,
    s.rx_name,
    s.control_iso_name,
    s.iso_name
"""


def gen_vp_dosage(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)
    records = cursor.fetchall()

    dosage_group = defaultdict(list)
    for rec in records:
        dosage = rec['dosage']
        dosage_group[dosage].append(rec)

    results = []
    for dosage, rx_list in dosage_group.items():
        results.append({
            'dosage': int(dosage),
            'num_ref_name': len(set(
                r['ref_name'] for r in rx_list
            )),
            'num_fold': sum([r['num_fold'] for r in rx_list])
        })

    results.sort(key=itemgetter('dosage'))

    dump_csv(DATA_FILE_PATH / 'table_vp_dosage.csv', results)
    dump_json(DATA_FILE_PATH / 'table_vp_dosage.json', results)
