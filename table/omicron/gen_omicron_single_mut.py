from statistics import median
from preset import DATA_FILE_PATH
from preset import dump_csv
from sql import row2dict
from preset import group_records_by

SQL = """
SELECT
    s.ref_name,
    rx.ab_name,
    s.fold_cmp,
    s.fold,
    s.control_potency,
    s.potency,
    mut.single_mut_name,
    mut.ref,
    mut.position,
    mut.amino_acid,
    mut.var_name
FROM
    susc_results_50_wt_view s,
    rx_mab_view rx,
    (
        select
            m.iso_name,
            m.single_mut_name,
            m.ref,
            m.position,
            m.amino_acid,
            group_concat(i.var_name, ', ') as var_name
        from
            isolate_mutations_single_s_mut_view m,
            isolate_mutations iso,
            isolates i
        where
            m.gene = iso.gene
            and
            m.position = iso.position
            and
            m.amino_acid = iso.amino_acid
            and
            iso.iso_name = i.iso_name
            and
            iso.iso_name in ('BA.1 Spike', 'BA.2 Spike', 'BA.1 Spike:+346K') and m.domain = 'RBD' and m.gene = 'S'
        group by
            m.iso_name,
            m.single_mut_name,
            m.ref,
            m.position,
            m.amino_acid
    ) as mut

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
ORDER BY
    mut.position
;
"""


def gen_omicron_single_mut(
        conn, folder=DATA_FILE_PATH / 'omicron'):

    cursor = conn.cursor()

    # for iso_name in ['BA.1 Spike', 'BA.2 Spike', 'BA.1 Spike:+346K']:

    #     sql = SQL.format(iso_name=iso_name)
    #     cursor.execute(sql)
    #     table = row2dict(cursor.fetchall())
    #     file_prefix = iso_name.split()[0]
    #     if iso_name == 'BA.1 Spike:+346K':
    #         file_prefix = 'BA.1.1'
    #     dump_csv(folder / f'{file_prefix}_single_mut.csv', table)

    cursor.execute(SQL)
    table = row2dict(cursor.fetchall())

    for i, mut_list in group_records_by(table, 'single_mut_name').items():
        for mab, mab_list in group_records_by(mut_list, 'ab_name').items():
            med = median([
                i['fold']
                for i in mab_list
            ])
            fold_values = [
                i['fold']
                for i in mab_list
            ]
            if med >=5:
                low_value = min([
                    i
                    for i in fold_values
                    if i >= 5
                ])
                high_value = max([
                    i
                    for i in fold_values
                    if i >= 5
                ])
                range = f'{low_value}-{high_value}'
            else:
                range = ''

            high_value = max(fold_values)
            low_value = min(fold_values)
            [
                r.update({
                    'median': med,
                    'range': range,
                })
                for r in mab_list
            ]


    dump_csv(folder / 'omicron_single_mut.csv', table)
