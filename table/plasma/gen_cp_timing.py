from preset import dump_csv
from preset import dump_json
from preset import DATA_FILE_PATH
from collections import defaultdict


SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    rx.timing,
    SUM(s.cumulative_count) num_fold
FROM
    susc_results_view s,
    rx_conv_plasma rx
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


def gen_cp_timing(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)
    records = cursor.fetchall()

    timing_group = defaultdict(list)
    for rec in records:
        timing = int(rec['timing'])
        if timing < 2:
            timing = '1 month'
        elif timing < 4:
            timing = '2-3 months'
        elif timing < 7:
            timing = '4-6 months'
        else:
            timing = '>6 months'
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

    dump_csv(DATA_FILE_PATH / 'table_cp_time.csv', timing_results)
    [
        i.update({
            'timing': i['timing'].replace('>', '&gt;')
        })
        for i in timing_results
    ]
    dump_json(DATA_FILE_PATH / 'table_cp_time.json', timing_results)
