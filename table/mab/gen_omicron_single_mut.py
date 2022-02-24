from preset import group_records_by
from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import row2dict

SQL = """
SELECT
    s.ref_name,
    rx.ab_name,
    s.fold_cmp,
    s.fold,
    control_potency,
    potency,
    mut.single_mut_name,
    mut.position
FROM
    susc_results_50_wt_view s,
    rx_mab_view rx,
    isolate_mutations_single_s_mut_view mut
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    rx.availability IS NOT NULL
    AND
    s.fold IS NOT NULL

    AND
    s.ref_name IN (
        SELECT
            DISTINCT susc.ref_name
        FROM
            susc_results susc,
            isolates iso,
            rx_mab_view mab
        WHERE
            susc.iso_name = iso.iso_name
            AND
            iso.var_name LIKE 'Omicron%'
            AND
            susc.ref_name = mab.ref_name
            AND
            susc.rx_name = mab.rx_name
            AND
            mab.availability IS NOT NULL
    )

    AND
    s.iso_name = mut.iso_name

    AND
    s.iso_name IN (
        SELECT
            iso_name
        FROM
            isolate_mutations_single_s_mut_view a
        WHERE EXISTS (
            SELECT 1
            FROM
                isolate_mutations b
            WHERE
                b.iso_name = '{variant_name}'
                AND
                a.position = b.position
                AND
                a.amino_acid = b.amino_acid
            )
            AND
            gene = 'S'
            AND
            domain = 'RBD'
            )
ORDER BY
    mut.position
;
"""


def gen_omicron_single_mut(
        conn, folder=DATA_FILE_PATH / 'mab'):

    cursor = conn.cursor()

    sql = SQL.format(variant_name='BA.2 Spike')
    cursor.execute(sql)
    table = row2dict(cursor.fetchall())
    dump_csv(folder / 'Omicron_BA_2_single_mut.csv', table)

    sql = SQL.format(variant_name='BA.1 Spike:+346K')
    cursor.execute(sql)
    table = row2dict(cursor.fetchall())
    dump_csv(folder / 'Omicron_BA_1_single_mut.csv', table)
