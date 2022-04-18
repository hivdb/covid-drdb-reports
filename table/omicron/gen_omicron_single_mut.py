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
    mut.ref,
    mut.position,
    mut.amino_acid
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
            iso.var_name IN ('Omicron/BA.1', 'Omicron/BA.2', 'Omicron/BA.1.1')
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
                isolate_mutations b,
                isolates c
            WHERE
                b.iso_name = c.iso_name
                AND
                c.var_name = '{var_name}'
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
        conn, folder=DATA_FILE_PATH / 'omicron'):

    cursor = conn.cursor()

    for var_name in ['Omicron/BA.1', 'Omicron/BA.2', 'Omicron/BA.1.1']:

        sql = SQL.format(var_name=var_name)
        cursor.execute(sql)
        table = row2dict(cursor.fetchall())
        file_prefix = var_name.split('/')[-1]
        dump_csv(folder / f'{file_prefix}_single_mut.csv', table)
