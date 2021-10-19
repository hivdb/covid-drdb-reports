from preset import DATA_FILE_PATH
from preset import dump_csv
from collections import defaultdict


CONTROL_50_SQL = """
SELECT
    wt.var_name control_iso_name,
    SUM(susc.cumulative_count) num_fold
FROM
    susc_results_50_wt_view susc,
    isolate_wildtype_view wt
WHERE
    susc.control_iso_name = wt.iso_name
GROUP BY
    susc.control_iso_name

UNION

SELECT
    CASE
        WHEN wt.var_name IS NOT NULL THEN
            wt.var_name
        ELSE
            wt.iso_name
    END control_iso_name,
    SUM(susc.cumulative_count) num_fold
FROM
    susc_results_50_view susc,
    isolate_non_wildtype_view wt
WHERE
    susc.control_iso_name = wt.iso_name
GROUP BY
    wt.var_name
"""

CONTROL_SQL = """
SELECT
    wt.var_name control_iso_name,
    SUM(susc.cumulative_count) num_fold
FROM
    susc_results_view susc,
    isolate_wildtype_view wt
WHERE
    susc.control_iso_name = wt.iso_name
GROUP BY
    susc.control_iso_name

UNION

SELECT
    CASE
        WHEN wt.var_name IS NOT NULL THEN
            wt.var_name
        ELSE
            wt.iso_name
    END control_iso_name,
    SUM(susc.cumulative_count) num_fold
FROM
    susc_results_view susc,
    isolate_non_wildtype_view wt
WHERE
    susc.control_iso_name = wt.iso_name
GROUP BY
    wt.var_name
"""


def gen_control(conn):
    save_path = DATA_FILE_PATH / 'fold' / 'control.csv'
    _gen_control(conn, save_path, CONTROL_SQL)

    save_path = DATA_FILE_PATH / 'fold' / 'control_50.csv'
    _gen_control(conn, save_path, CONTROL_50_SQL)


def _gen_control(conn, save_path, sql):
    cursor = conn.cursor()
    cursor.execute(sql)

    records = cursor.fetchall()

    control_name_group = defaultdict(list)
    for rec in records:
        control = rec['control_iso_name']
        control_name_group[control].append(rec)

    control_name_results = []
    for control, rx_list in control_name_group.items():
        control_name_results.append({
            'control': control,
            'num_fold': sum([r['num_fold'] for r in rx_list])
        })

    dump_csv(save_path, control_name_results)
