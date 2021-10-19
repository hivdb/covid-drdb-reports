from preset import DATA_FILE_PATH
from preset import dump_csv

NULL_FOLD_SQL = """
SELECT
    susc.ref_name,
    susc.rx_name,
    susc.resistance_level,
    susc.potency_type,
    SUM(susc.cumulative_count) num_fold
FROM
    susc_results_view susc
WHERE
    fold IS NULL
GROUP BY
    susc.ref_name,
    susc.rx_name,
    susc.resistance_level,
    susc.potency_type
    ;
"""


def gen_null_fold(
        conn,
        save_path=DATA_FILE_PATH / 'fold' / 'null_fold.csv'):

    cursor = conn.cursor()

    cursor.execute(NULL_FOLD_SQL)

    results = []
    for rec in cursor.fetchall():
        results.append({
            'ref_name': rec['ref_name'],
            'rx_name': rec['rx_name'],
            'resistance': rec['resistance_level'],
            'potency': rec['potency_type'],
            'num_fold': rec['num_fold']
        })

    dump_csv(save_path, results)
