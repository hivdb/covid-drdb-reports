from preset import dump_csv
from preset import DATA_FILE_PATH


SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.iso_name,
    var.as_wildtype,
    rx.timing,
    rx.severity,
    rx.infected_var_name,
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
    rx.infected_var_name = iso.var_name
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


def gen_cp_summary(conn):
    cursor = conn.cursor()
    sql = SQL.format(
        susc_results_type='susc_results_indiv_view'
    )

    cursor.execute(sql)
    records = cursor.fetchall()

    num_fold_results = sum([r['num_fold'] for r in records])
    num_ref_name = len(set([r['ref_name'] for r in records]))

    sql = SQL.format(
        susc_results_type='susc_results_aggr_view'
    )

    cursor.execute(sql)
    aggre_records = cursor.fetchall()

    aggre_num_fold_results = sum([r['num_fold'] for r in aggre_records])
    aggre_num_ref_study = len(set([r['ref_name'] for r in aggre_records]))

    all_num_ref_name = num_fold_results + aggre_num_fold_results
    all_num_fold = len(set([r['ref_name'] for r in records + aggre_records]))

    result = [{
        'indiv_num_fold': num_fold_results,
        'indiv_num_fold_ref_name': num_ref_name,
        'aggre_num_fold': aggre_num_fold_results,
        'aggre_num_fold_ref_name': aggre_num_ref_study,
        'all_num_fold': all_num_fold,
        'all_num_ref_name': all_num_ref_name,
    }]
    save_path = DATA_FILE_PATH / 'cp' / 'summary_cp.csv'
    dump_csv(save_path, result)
