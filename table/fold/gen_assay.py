from preset import DATA_FILE_PATH
from preset import dump_csv
from collections import defaultdict


ASSAY_SQL = """
SELECT
    susc.assay_name,
    SUM(susc.cumulative_count) AS num_fold
FROM
    susc_results susc,
    isolate_wildtype_view wt
WHERE
    susc.control_iso_name = wt.iso_name
GROUP BY
    susc.assay_name
    ;
"""


ASSAY_50_SQL = """
SELECT
    susc.assay_name,
    SUM(susc.cumulative_count) AS num_fold
FROM
    susc_results_50_wt_view AS susc
GROUP BY
    susc.assay_name
    ;
"""


def gen_assay(conn):
    sql = ASSAY_50_SQL
    save_path = DATA_FILE_PATH / 'fold' / 'assay_50.csv'
    _gen_assay(conn, save_path, sql)

    sql = ASSAY_SQL
    save_path = DATA_FILE_PATH / 'fold' / 'assay.csv'
    _gen_assay(conn, save_path, sql)


def _gen_assay(conn, save_path, sql):
    cursor = conn.cursor()
    cursor.execute(sql)

    records = cursor.fetchall()

    assay_group = defaultdict(list)
    for rec in records:
        assay = rec['assay_name']
        assay_group[assay].append(rec)

    assay_results = []
    for assay, rx_list in assay_group.items():
        assay_results.append({
            'assay': assay,
            'num_fold': sum([r['num_fold'] for r in rx_list])
        })

    dump_csv(save_path, assay_results)
