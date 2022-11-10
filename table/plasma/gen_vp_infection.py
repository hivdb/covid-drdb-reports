from preset import dump_csv
from preset import dump_json
from preset import DATA_FILE_PATH
from collections import defaultdict
from operator import itemgetter

SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    rx.as_wildtype,
    rx.var_name infection,
    SUM(s.cumulative_count) num_fold
FROM
    susc_results_view s,
    rx_vacc_plasma_infect_var_view rx
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


def gen_vp_infection(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)
    records = cursor.fetchall()

    infection_group = defaultdict(list)
    for rec in records:
        infection = rec['infection']
        as_wildtype = rec['as_wildtype']
        if as_wildtype:
            infection = 'Wiletype'

        if infection:
            infection_group['infected'].append(rec)
            infection_group[infection].append(rec)
        else:
            infection_group['uninfected'].append(rec)

    infection_results = []
    for infection, rx_list in infection_group.items():
        infection_results.append({
            'infection': infection,
            'num_ref_name': len(set(
                r['ref_name'] for r in rx_list
            )),
            'num_fold': sum([r['num_fold'] for r in rx_list])
        })
    infection_results.sort(key=itemgetter('infection'))
    dump_csv(DATA_FILE_PATH / 'table_vp_infection.csv', infection_results)
    dump_json(DATA_FILE_PATH / 'table_vp_infection.json', infection_results)
