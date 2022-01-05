from preset import DATA_FILE_PATH
from preset import dump_csv
import re


SUMMARY_SQL = """
SELECT DISTINCT
    s.ref_name,
    iso.pattern
FROM
    susc_results_view s,
    {rx_type} rx,
    {iso_type} iso
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.iso_name = iso.iso_name
    AND
    iso.var_name = 'Omicron'
;
"""


def gen_omicron_s_pattern(
        conn,
        csv_save_path=DATA_FILE_PATH / 'mab' / 'omicron_s_pattern.csv'):
    cursor = conn.cursor()

    sql = SUMMARY_SQL.format(
        rx_type='rx_mab_view',
        iso_type='isolate_mutations_combo_s_mut_view',
    )

    cursor.execute(sql)

    results = []
    mut_pos_list = set()
    pos_ref_aa_dict = {}
    for row in cursor.fetchall():
        ref_name = row['ref_name']
        pattern = row['pattern']
        results.append([ref_name, pattern])

        mut_list = pattern.split('+')
        for mut in mut_list:
            ref_aa, pos, mut_aa = re.match(
                r'^(\D+)([\d\-]+)(\D+)$', mut).groups()

            if '-' in pos:
                pos = pos.split('-', 1)[0]
            pos = int(pos)

            mut_pos_list.add(pos)
            pos_ref_aa_dict[pos] = ref_aa

    mut_pos_list = sorted(list(mut_pos_list))

    headers = ['ref_name']
    for pos in mut_pos_list:
        ref_aa = pos_ref_aa_dict[pos]
        headers.append('{}{}'.format(ref_aa, pos))

    summary_results = []
    for ref_name, pattern in results:
        result = {
            'ref_name': ref_name
        }

        mut_list = pattern.split('+')
        for mut in mut_list:
            ref_aa, pos, mut_aa = re.match(
                r'^(\D+)([\d\-]+)(\D+)$', mut).groups()
            if '-' in pos:
                pos = pos.split('-', 1)[0]

            result['{}{}'.format(ref_aa, pos)] = mut_aa

        summary_results.append(result)

    dump_csv(csv_save_path, summary_results, headers=headers)
