from preset import dump_csv
from preset import DATA_FILE_PATH
from preset import round_number
from collections import defaultdict
from operator import itemgetter
from statistics import median
from resistancy import is_susc
from resistancy import is_partial_resistant
from resistancy import is_resistant

SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.fold,
    rx.vaccine_name,
    SUM(s.cumulative_count) num_fold
FROM
    susc_results as s,
    rx_vacc_plasma as rx
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


def gen_vp_vaccine(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)
    records = cursor.fetchall()
    num_records = sum([r['num_fold'] for r in records])

    vaccine_group = defaultdict(list)
    for rec in records:
        vaccine = rec['vaccine_name']
        vaccine_group[vaccine].append(rec)

    vaccine_results = []
    for vaccine, rx_list in vaccine_group.items():

        all_fold = [[r['fold']] * r['num_fold'] for r in rx_list]
        all_fold = [r for j in all_fold for r in j]
        num_fold = len(all_fold)

        s_fold = [r for r in all_fold if is_susc(r)]
        i_fold = [r for r in all_fold if is_partial_resistant(r)]
        r_fold = [r for r in all_fold if is_resistant(r)]

        num_s_fold = len(s_fold)
        num_i_fold = len(i_fold)
        num_r_fold = len(r_fold)
        median_fold = median(all_fold)

        vaccine_results.append({
            'vaccine': vaccine,
            'num_ref_name': len(set(
                r['ref_name'] for r in rx_list
            )),
            'num_fold': num_fold,
            'S': num_s_fold,
            'I': num_i_fold,
            'R': num_r_fold,
            'median_fold': median_fold,
            'percent': round_number(num_fold / num_records * 100),
        })
    vaccine_results.sort(key=itemgetter('vaccine'))
    save_path = DATA_FILE_PATH / 'vp' / 'summary_vp_vaccine.csv'
    dump_csv(save_path, vaccine_results)
