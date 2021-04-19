from preset import dump_csv
from preset import DATA_FILE_PATH
from variant.preset import CONTROL_VARIANTS_SQL


VP_SQL = """
SELECT
    rx.vaccine_name as vaccine_name,
    COUNT(s.cumulative_count) as samples
FROM
    susc_results as s,
    rx_immu_plasma as rx
ON
    s.ref_name = rx.ref_name
    AND s.rx_name = rx.rx_name
WHERE
    s.control_variant_name in {control_variants}
    AND s.fold IS NOT NULL
GROUP BY rx.vaccine_name
""".format(control_variants=CONTROL_VARIANTS_SQL)


def gen_table_vp(
        conn, csv_save_path=DATA_FILE_PATH / "summary_vp.csv"):

    cursor = conn.cursor()
    cursor.execute(VP_SQL)

    results = []
    for item in cursor.fetchall():
        results.append({
            'Vaccine': item['vaccine_name'],
            'Samples': item['samples']
        })

    dump_csv(csv_save_path, results)
