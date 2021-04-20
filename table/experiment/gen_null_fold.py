from preset import DATA_FILE_PATH
from preset import dump_csv

SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.resistance_level,
    SUM(s.cumulative_count) as sample
FROM
    susc_results AS s
WHERE
    fold IS NULL
GROUP BY
    s.ref_name,
    s.rx_name,
    s.resistance_level
    ;
"""


def gen_null_fold(conn, save_path=DATA_FILE_PATH / 'summary_null_fold.csv'):

    cursor = conn.cursor()

    cursor.execute(SQL)

    results = []
    for rec in cursor.fetchall():
        results.append({
            'ref_name': rec['ref_name'],
            'rx_name': rec['rx_name'],
            'resistance': rec['resistance_level'],
            'sample': rec['sample']
        })

    dump_csv(save_path, results)
