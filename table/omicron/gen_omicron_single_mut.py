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
    mut.domain = 'RBD'
    AND
    EXISTS (
        SELECT
        1
        FROM
            isolate_mutations b
        WHERE
            b.iso_name = '{iso_name}'
            AND
            mut.gene = b.gene
            AND
            mut.position = b.position
            AND
            mut.amino_acid = b.amino_acid
            AND
            b.gene = 'S'
        )
ORDER BY
    mut.position
;
"""


def gen_omicron_single_mut(
        conn, folder=DATA_FILE_PATH / 'omicron'):

    cursor = conn.cursor()

    for iso_name in ['BA.1 Spike', 'BA.2 Spike', 'BA.1 Spike:+346K']:

        sql = SQL.format(iso_name=iso_name)
        cursor.execute(sql)
        table = row2dict(cursor.fetchall())
        file_prefix = iso_name.split()[0]
        if iso_name == 'BA.1 Spike:+346K':
            file_prefix = 'BA.1.1'
        dump_csv(folder / f'{file_prefix}_single_mut.csv', table)
