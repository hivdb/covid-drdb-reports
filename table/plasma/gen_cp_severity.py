from preset import dump_csv
from preset import DATA_FILE_PATH
from collections import defaultdict


SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    rx.severity,
    SUM(s.cumulative_count) as num_fold
FROM
    susc_results_view s,
    rx_conv_plasma rx
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


def gen_cp_severity(conn):
    cursor = conn.cursor()
    cursor.execute(SQL)
    records = cursor.fetchall()

    severity_group = defaultdict(list)
    for rec in records:
        severity = rec['severity']
        severity_group[severity].append(rec)

    severiy_results = []
    for severity, rx_list in severity_group.items():
        severiy_results.append({
            'severity': severity,
            'num_ref_name': len(set(
                r['ref_name'] for r in rx_list
            )),
            'num_fold': sum([r['num_fold'] for r in rx_list])
        })
    save_path = DATA_FILE_PATH / 'cp' / 'summary_cp_severity.csv'
    dump_csv(save_path, severiy_results)
