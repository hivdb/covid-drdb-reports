from preset import DATA_FILE_PATH
from preset import dump_csv
from variant.preset import filter_by_variant
from collections import defaultdict


SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.control_iso_name,
    s.iso_name,
    s.assay,
    SUM(s.cumulative_count) AS samples
FROM
    susc_results AS s
WHERE
    s.inhibition_pcnt != 90
GROUP BY
    s.ref_name,
    s.rx_name,
    s.control_iso_name,
    s.iso_name
    ;
"""


def gen_exp(conn):
    cursor = conn.cursor()
    cursor.execute(SQL)

    records = cursor.fetchall()
    records = filter_by_variant(records)

    assay_group = defaultdict(list)
    for rec in records:
        assay = rec['assay']
        assay_group[assay].append(rec)

    assay_results = []
    for assay, rx_list in assay_group.items():
        assay_results.append({
            'Assay': assay,
            'Samples': sum([r['samples'] for r in rx_list])
        })

    save_path = DATA_FILE_PATH / 'summary_exp_assay.csv'
    dump_csv(save_path, assay_results)

    control_name_group = defaultdict(list)
    for rec in records:
        control = rec['control_iso_name']
        control_name_group[control].append(rec)

    control_name_results = []
    for control, rx_list in control_name_group.items():
        control_name_results.append({
            'Control': control,
            'Samples': sum([r['samples'] for r in rx_list])
        })

    save_path = DATA_FILE_PATH / 'summary_exp_control.csv'
    dump_csv(save_path, control_name_results)
