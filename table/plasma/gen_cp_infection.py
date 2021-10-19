from preset import dump_csv
from preset import DATA_FILE_PATH
from collections import defaultdict


SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    rx.as_wildtype,
    rx.var_name infection,
    SUM(s.cumulative_count) num_fold
FROM
    susc_results_view s,
    rx_conv_plasma_infect_var_view rx
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


def gen_cp_infection(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)
    records = cursor.fetchall()

    infection_group = defaultdict(list)
    for rec in records:
        infection = rec['infection']
        as_wildtype = rec['as_wildtype']
        if as_wildtype:
            infection = 'wt'

        if infection:
            infection_group['infected'].append(rec)
            infection_group[infection].append(rec)
        else:
            infection_group['unknown'].append(rec)

    infection_results = []
    for infection, rx_list in infection_group.items():
        infection_results.append({
            'infection': infection,
            'num_ref_name': len(set(
                r['ref_name'] for r in rx_list
            )),
            'num_fold': sum([r['num_fold'] for r in rx_list])
        })
    save_path = DATA_FILE_PATH / 'cp' / 'summary_cp_infection.csv'
    dump_csv(save_path, infection_results)
