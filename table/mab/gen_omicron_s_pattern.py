from preset import DATA_FILE_PATH
from preset import dump_csv
import re
from collections import defaultdict


SUMMARY_SQL = """
SELECT DISTINCT
    s.ref_name,
    iso.pattern,
    CASE
        WHEN iso_meta.gisaid_id IS NOT NULL THEN
            iso_meta.iso_name
        WHEN iso_meta.genbank_accn IS NOT NULL THEN
            iso_meta.iso_name
        ELSE
            ''
    END iso_name
FROM
    susc_results_view s,
    {rx_type} rx,
    {iso_type} iso,
    isolates iso_meta
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.iso_name = iso.iso_name
    AND
    iso.var_name = 'Omicron'
    AND
    iso.iso_name = iso_meta.iso_name
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
    pos_mut_stat = {}
    for row in cursor.fetchall():
        ref_name = row['ref_name']
        iso_name = row['iso_name']
        pattern = row['pattern']
        results.append([ref_name, iso_name, pattern])

        mut_list = pattern.split('+')
        for mut in mut_list:
            ref_aa, pos, mut_aa = re.match(
                r'^(\D+)([\d\-]+)(\D+)$', mut).groups()

            if '-' in pos:
                pos = pos.split('-', 1)[0]
            pos = int(pos)

            mut_pos_list.add(pos)
            pos_ref_aa_dict[pos] = ref_aa

            if pos not in pos_mut_stat:
                pos_mut_stat[pos] = defaultdict(int)
            pos_mut_stat[pos][mut_aa] += 1

    mut_pos_list = sorted(list(mut_pos_list))

    pos_most_common_mut_aa = {}
    for pos, mut_aa_stat in pos_mut_stat.items():
        mut_aa_stat = [(k, v) for k, v in mut_aa_stat.items() if k]
        most_common_mut_aa = sorted(mut_aa_stat, key=lambda x: x[1])[-1][0]
        pos_most_common_mut_aa[pos] = most_common_mut_aa

    headers = ['ref_name', 'iso_name']
    for pos in mut_pos_list:
        ref_aa = pos_ref_aa_dict[pos]
        mut_aa = pos_most_common_mut_aa[pos]
        headers.append('{}{}{}'.format(ref_aa, pos, mut_aa))

    summary_results = []
    for ref_name, iso_name, pattern in results:
        result = {
            'ref_name': ref_name,
            'iso_name': iso_name
        }

        mut_list = pattern.split('+')
        for mut in mut_list:
            ref_aa, pos, mut_aa = re.match(
                r'^(\D+)([\d\-]+)(\D+)$', mut).groups()
            if '-' in pos:
                pos = pos.split('-', 1)[0]

            pos = int(pos)

            most_common_mut_aa = pos_most_common_mut_aa[pos]

            result['{}{}{}'.format(ref_aa, pos, most_common_mut_aa)] = (
                '-' if most_common_mut_aa == mut_aa else mut_aa
            )

        summary_results.append(result)

    dump_csv(csv_save_path, summary_results, headers=headers)
