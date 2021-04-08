from preset import dump_csv
from preset import DATA_FILE_PATH


VP_SQL = """
SELECT
    b.vaccine_name as vaccine_name,
    count(a.cumulative_count) as samples
FROM
    susc_results as a,
    rx_immu_plasma as b
ON
    a.ref_name = b.ref_name AND a.rx_name = b.rx_name
GROUP BY b.vaccine_name
"""


def gen_table_vp(
        conn, csv_save_path=DATA_FILE_PATH / "table_vp_figure.csv"):

    cursor = conn.cursor()
    cursor.execute(VP_SQL)

    results = []
    for item in cursor.fetchall():
        results.append({
            'Vaccine': item['vaccine_name'],
            'Samples': item['samples']
        })

    dump_csv(csv_save_path, results)
