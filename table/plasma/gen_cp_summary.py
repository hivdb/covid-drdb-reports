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
    var.as_wildtype,
    rx.timing,
    rx.severity,
    rx.infected_iso_name,
    iso.var_name AS infection,
    SUM(s.cumulative_count) as num_results
FROM
    ({susc_results}) as s,
    rx_conv_plasma as rx,
    isolates as iso,
    variants as var
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.iso_name = iso.iso_name
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
# -- s.control_iso_name in ({control_variants})
# -- AND


def gen_cp_summary(conn):
    cursor = conn.cursor()
    sql = SQL.format(
        susc_results=INDIVIDUAL_SAMPLE_SQL,
        control_variants=CONTROL_VARIANTS_SQL
    )

    records = cursor.execute(sql)
    records = filter_by_variant(records)

    num_fold_results = sum([r['num_results'] for r in records])
    num_study = len(set([r['ref_name'] for r in records]))

    sql = SQL.format(
        susc_results=AGGREGATED_SUSC_VIEW_SQL,
        control_variants=CONTROL_VARIANTS_SQL
    )

    aggre_records = cursor.execute(sql)
    aggre_records = filter_by_variant(aggre_records)

    aggre_num_fold_results = sum([r['num_results'] for r in aggre_records])
    aggre_num_study = len(set([r['ref_name'] for r in aggre_records]))

    result = [{
        'indiv_num_results': num_fold_results,
        'indiv_num_results_study': num_study,
        'aggre_num_results': aggre_num_fold_results,
        'aggre_num_results_study': aggre_num_study,
    }]
    save_path = DATA_FILE_PATH / 'summary_cp.csv'
    dump_csv(save_path, result)

    records += aggre_records

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
            'Timing': timing,
            'References': len(set(
                r['ref_name'] for r in rx_list
            )),
            'Results': sum([r['num_results'] for r in rx_list])
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
            'References': len(set(
                r['ref_name'] for r in rx_list
            )),
            'Results': sum([r['num_results'] for r in rx_list])
        })
    save_path = DATA_FILE_PATH / 'summary_cp_severity.csv'
    dump_csv(save_path, severiy_results)

    infection_group = defaultdict(list)
    for rec in records:
        infection = rec['infection']
        as_wildtype = rec['as_wildtype']
        if as_wildtype:
            infection = 'wt'
        else:
            infection = infection.split()[0]
            infection = infection.split('/')[0]
        infection_group[infection].append(rec)

    infection_results = []
    for infection, rx_list in infection_group.items():
        infection_results.append({
            'Infection': infection,
            'References': len(set(
                r['ref_name'] for r in rx_list
            )),
            'Results': sum([r['num_results'] for r in rx_list])
        })
    save_path = DATA_FILE_PATH / 'summary_cp_infection.csv'
    dump_csv(save_path, infection_results)
