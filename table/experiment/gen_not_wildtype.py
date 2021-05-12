from variant.preset import CONTROL_VARIANTS_SQL
from preset import DATA_FILE_PATH
from preset import dump_csv

SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.control_iso_name,
    s.iso_name,
    SUM(s.cumulative_count) AS sample
FROM
    susc_results AS s
WHERE
    s.inhibition_pcnt != 90
    AND
    s.fold IS NOT NULL
    AND
    control_iso_name NOT IN {control_variants}
GROUP BY
    s.ref_name,
    s.rx_name,
    s.control_iso_name,
    s.iso_name
    ;
""".format(control_variants=CONTROL_VARIANTS_SQL)


def gen_not_wildtype(
        conn,
        save_path=DATA_FILE_PATH / 'summary_not_wildtype.csv'):

    cursor = conn.cursor()
    cursor.execute(SQL)

    results = []
    for rec in cursor.fetchall():
        results.append({
            'ref_name': rec['ref_name'],
            'rx_name': rec['rx_name'],
            'control_variant': rec['control_iso_name'],
            'variant': rec['iso_name'],
            'sample': rec['sample']
        })

    dump_csv(save_path, results)
