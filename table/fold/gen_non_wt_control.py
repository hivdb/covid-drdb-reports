from preset import DATA_FILE_PATH
from preset import dump_csv

NON_WT_CONTROL_SQL = """
SELECT
    susc.ref_name,
    susc.rx_name,
    CASE
        WHEN wt.var_name IS NOT NULL THEN
            wt.var_name
        ELSE
            wt.iso_name
    END control_iso_name,
    susc.iso_name,
    susc.potency_type,
    SUM(susc.cumulative_count) AS num_fold
FROM
    susc_results_view susc,
    isolate_non_wildtype_view wt
WHERE
    susc.control_iso_name = wt.iso_name
GROUP BY
    susc.ref_name,
    susc.rx_name,
    susc.control_iso_name,
    susc.iso_name,
    susc.potency_type
;
"""


def gen_non_wt_control(
        conn,
        save_path=DATA_FILE_PATH / 'fold' / 'non_wildtype.csv'):

    cursor = conn.cursor()
    cursor.execute(NON_WT_CONTROL_SQL)

    results = []
    for rec in cursor.fetchall():
        results.append({
            'ref_name': rec['ref_name'],
            'rx_name': rec['rx_name'],
            'control_var': rec['control_iso_name'],
            'iso_name': rec['iso_name'],
            'potency_type': rec['potency_type'],
            'num_fold': rec['num_fold']
        })

    dump_csv(save_path, results)
