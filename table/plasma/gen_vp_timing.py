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
    rx.timing,
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


def gen_vp_timing(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)
    records = cursor.fetchall()

    timing_group = defaultdict(list)
    for rec in records:
        timing = int(rec['timing'])
        if timing < 2:
            timing = '1'
        elif timing < 4:
            timing = '2-3'
        elif timing < 7:
            timing = '4-6'
        else:
            timing = '>6'
        timing_group[timing].append(rec)

    timing_results = []
    for timing, rx_list in timing_group.items():
        timing_results.append({
            'timing': timing,
            'num_ref_name': len(set(
                r['ref_name'] for r in rx_list
            )),
            'num_fold': sum([r['num_fold'] for r in rx_list])
        })

    timing_results.sort(key=lambda x: x['timing'] if x['timing'] else 0)

    save_path = DATA_FILE_PATH / 'vp' / 'summary_vp_timing.csv'
    dump_csv(save_path, timing_results)
