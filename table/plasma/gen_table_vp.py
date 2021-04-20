from preset import dump_csv
from preset import DATA_FILE_PATH
from collections import defaultdict
from operator import itemgetter
from variant.preset import CONTROL_VARIANTS_SQL
from variant.preset import filter_by_variant


VP_SQL = """
SELECT
    s.variant_name,
    rx.vaccine_name,
    SUM(s.cumulative_count) as samples
FROM
    susc_results as s,
    rx_immu_plasma as rx
ON
    s.ref_name = rx.ref_name
    AND s.rx_name = rx.rx_name
WHERE
    s.control_variant_name in {control_variants}
    AND s.fold IS NOT NULL
GROUP BY
    s.variant_name,
    rx.vaccine_name
""".format(control_variants=CONTROL_VARIANTS_SQL)


def gen_table_vp(
        conn, csv_save_path=DATA_FILE_PATH / "summary_vp.csv"):

    cursor = conn.cursor()
    cursor.execute(VP_SQL)

    results = []
    records = cursor.fetchall()
    records = filter_by_variant(records)

    vaccine_groups = defaultdict(list)
    for rec in records:
        vaccine = rec['vaccine_name']
        vaccine_groups[vaccine].append(rec)

    for vaccine, rx_list in vaccine_groups.items():
        results.append({
            'Vaccine': vaccine,
            'Samples': sum([r['samples'] for r in rx_list])
        })

    results.sort(key=itemgetter('Vaccine'))

    dump_csv(csv_save_path, results)
