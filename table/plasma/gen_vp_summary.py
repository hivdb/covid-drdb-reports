from preset import dump_csv
from preset import DATA_FILE_PATH
from preset import round_number
from collections import defaultdict
from statistics import median
from resistancy import is_susc
from resistancy import is_partial_resistant
from resistancy import is_resistant

SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.iso_name,
    s.fold,
    rx.timing,
    rx.dosage,
    rx.vaccine_name,
    iso.var_name,
    sub.subject_species species,
    SUM(s.cumulative_count) num_fold
FROM
    {susc_results_type} s,
    rx_vacc_plasma rx,
    subjects sub,
    isolates iso
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    rx.ref_name = sub.ref_name
    AND
    rx.subject_name = sub.subject_name
    AND
    s.iso_name = iso.iso_name
    AND
    s.fold IS NOT NULL
GROUP BY
    s.ref_name,
    s.rx_name,
    s.control_iso_name,
    s.iso_name
"""


def gen_vp_summary(conn):
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
    save_path = DATA_FILE_PATH / 'vp' / 'summary_vp.csv'
    dump_csv(save_path, result)
