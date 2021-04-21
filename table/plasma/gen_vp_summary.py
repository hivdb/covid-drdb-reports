from preset import dump_csv
from preset import DATA_FILE_PATH
from preset import round_number
from collections import defaultdict
from statistics import median
from operator import itemgetter
from variant.preset import CONTROL_VARIANTS_SQL
from variant.preset import filter_by_variant
from lookup_view import INDIVIDUAL_SAMPLE_SQL
from lookup_view import AGGREGATED_SUSC_VIEW_SQL
from resistancy import is_susc
from resistancy import is_partial_resistant
from resistancy import is_resistant

SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.variant_name,
    s.fold,
    rx.timing,
    rx.dosage,
    rx.vaccine_name,
    SUM(s.cumulative_count) as samples
FROM
    ({susc_results}) as s,
    rx_vacc_plasma as rx
ON
    s.ref_name = rx.ref_name
    AND s.rx_name = rx.rx_name
WHERE
    s.fold IS NOT NULL
GROUP BY
    s.ref_name,
    s.rx_name,
    s.control_variant_name,
    s.variant_name
"""
# -- s.control_variant_name in {control_variants}
# -- AND


def gen_vp_summary(conn):
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
    save_path = DATA_FILE_PATH / 'summary_vp.csv'
    dump_csv(save_path, result)

    records += aggre_records
    num_records = sum([r['samples'] for r in records])

    timing_dosage_group = defaultdict(list)
    for rec in records:
        timing = rec['timing']
        dosage = rec['dosage']
        timing_dosage_group[(timing, dosage)].append(rec)

    timing_results = []
    for (timing, dosage), rx_list in timing_dosage_group.items():
        timing_results.append({
            'Timing': timing,
            'Dosage': dosage,
            'Samples': sum([r['samples'] for r in rx_list])
        })

    save_path = DATA_FILE_PATH / 'summary_vp_timing.csv'
    dump_csv(save_path, timing_results)

    vaccine_group = defaultdict(list)
    for rec in records:
        vaccine = rec['vaccine_name']
        vaccine_group[vaccine].append(rec)

    vaccine_results = []
    for vaccine, rx_list in vaccine_group.items():

        all_fold = [[r['fold']] * r['samples'] for r in rx_list]
        all_fold = [r for j in all_fold for r in j]
        samples = len(all_fold)

        s_fold = [r for r in all_fold if is_susc(r)]
        i_fold = [r for r in all_fold if is_partial_resistant(r)]
        r_fold = [r for r in all_fold if is_resistant(r)]

        num_s_fold = len(s_fold)
        num_i_fold = len(i_fold)
        num_r_fold = len(r_fold)
        median_fold = median(all_fold)

        vaccine_results.append({
            'Vaccine': vaccine,
            'Samples': samples,
            'S': num_s_fold,
            'I': num_i_fold,
            'R': num_r_fold,
            'median fold': median_fold,
            'Percent': round_number(samples / num_records * 100),
        })
    save_path = DATA_FILE_PATH / 'summary_vp_vaccine.csv'
    dump_csv(save_path, vaccine_results)
