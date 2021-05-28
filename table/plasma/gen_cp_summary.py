from preset import dump_csv
from preset import DATA_FILE_PATH
from collections import defaultdict
from operator import itemgetter
from variant.preset import CONTROL_VARIANTS_SQL
from variant.preset import filter_by_variant
from susceptibility import INDIVIDUAL_SAMPLE_SQL
from susceptibility import AGGREGATED_SUSC_VIEW_SQL


SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.iso_name,
    rx.timing,
    rx.severity,
    rx.infected_iso_name as infection,
    SUM(s.cumulative_count) as samples
FROM
    ({susc_results}) as s,
    rx_conv_plasma as rx
ON
    s.ref_name = rx.ref_name
    AND s.rx_name = rx.rx_name
WHERE
    s.fold IS NOT NULL
GROUP BY
    s.ref_name,
    s.rx_name,
    s.control_iso_name,
    s.iso_name
"""
# -- s.control_iso_name in {control_variants}
# -- AND


def gen_cp_summary(conn):
    cursor = conn.cursor()
    sql = SQL.format(
        susc_results=INDIVIDUAL_SAMPLE_SQL,
        control_variants=CONTROL_VARIANTS_SQL
    )

    records = cursor.execute(sql)
    records = filter_by_variant(records)

    num_exp = sum([r['samples'] for r in records])
    num_study = len(set([r['ref_name'] for r in records]))

    sql = SQL.format(
        susc_results=AGGREGATED_SUSC_VIEW_SQL,
        control_variants=CONTROL_VARIANTS_SQL
    )

    aggre_records = cursor.execute(sql)
    aggre_records = filter_by_variant(aggre_records)

    aggre_num_exp = sum([r['samples'] for r in aggre_records])
    aggre_num_study = len(set([r['ref_name'] for r in aggre_records]))

    result = [{
        'indiv_samples': num_exp,
        'indiv_samples_study': num_study,
        'aggre_samples': aggre_num_exp,
        'aggre_samples_study': aggre_num_study,
    }]
    save_path = DATA_FILE_PATH / 'summary_cp.csv'
    dump_csv(save_path, result)

    records += aggre_records

    timing_group = defaultdict(list)
    for rec in records:
        timing = rec['timing']
        timing_group[timing].append(rec)

    timing_results = []
    for timing, rx_list in timing_group.items():
        timing_results.append({
            'Timing': timing,
            'Samples': sum([r['samples'] for r in rx_list])
        })

    timing_results.sort(key=lambda x: x['Timing'] if x['Timing'] else 0)

    save_path = DATA_FILE_PATH / 'summary_cp_timing.csv'
    dump_csv(save_path, timing_results)

    severity_group = defaultdict(list)
    for rec in records:
        severity = rec['severity']
        severity_group[severity].append(rec)

    severiy_results = []
    for severity, rx_list in severity_group.items():
        severiy_results.append({
            'Severity': severity,
            'Samples': sum([r['samples'] for r in rx_list])
        })
    save_path = DATA_FILE_PATH / 'summary_cp_severity.csv'
    dump_csv(save_path, severiy_results)

    infection_group = defaultdict(list)
    for rec in records:
        infection = rec['infection']
        infection_group[infection].append(rec)

    infection_results = []
    for infection, rx_list in infection_group.items():
        infection_results.append({
            'Infection': infection,
            'Samples': sum([r['samples'] for r in rx_list])
        })
    save_path = DATA_FILE_PATH / 'summary_cp_infection.csv'
    dump_csv(save_path, infection_results)
