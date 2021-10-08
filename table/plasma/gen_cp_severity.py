from preset import dump_csv
from preset import DATA_FILE_PATH
from collections import defaultdict


SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.iso_name,
    var.as_wildtype,
    rx.timing,
    rx.severity,
    rx.infected_iso_name,
    iso.var_name AS infection,
    SUM(s.cumulative_count) as num_fold
FROM
    {susc_results_type} as s,
    rx_conv_plasma as rx,
    isolates as iso,
    variants as var
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.iso_name = iso.iso_name
    AND
    s.fold IS NOT NULL
    AND
    iso.var_name = var.var_name
GROUP BY
    s.ref_name,
    s.rx_name,
    s.control_iso_name,
    s.iso_name
"""


def gen_cp_severity(conn):
    cursor = conn.cursor()
    sql = SQL.format(
        susc_results_type='susc_results_indiv_view'
    )

    cursor.execute(sql)
    records = cursor.fetchall()

    sql = SQL.format(
        susc_results_type='susc_results_aggr_view'
    )

    cursor.execute(sql)
    aggre_records = cursor.fetchall()

    records += aggre_records

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
