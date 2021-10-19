from preset import dump_csv
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
            'dosage': dosage,
            'num_ref_name': len(set(
                r['ref_name'] for r in rx_list
            )),
            'num_fold': sum([r['num_fold'] for r in rx_list])
        })
    save_path = DATA_FILE_PATH / 'vp' / 'summary_vp_dosage.csv'
    dump_csv(save_path, results)
