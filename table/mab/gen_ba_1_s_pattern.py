from enum import unique
from preset import DATA_FILE_PATH
from preset import dump_csv
import re
from collections import defaultdict


SUMMARY_SQL = """
SELECT DISTINCT
    s.ref_name,
    ctl.var_name control_name,
    test.pattern,
    CASE
        WHEN test_meta.gisaid_id IS NOT NULL THEN
            test_meta.iso_name
        WHEN test_meta.genbank_accn IS NOT NULL THEN
            test_meta.iso_name
        ELSE
            ''
    END iso_name
FROM
    susc_results_view s,
    {rx_type} rx,
    isolates ctl,
    {iso_type} test,
    isolates test_meta
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name

    AND
    s.control_iso_name = ctl.iso_name

    AND
    s.iso_name = test.iso_name
    AND
    test.var_name = 'Omicron/BA.1'
    AND
    test.iso_name = test_meta.iso_name
;
"""


def gen_ba_1_s_pattern(
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
        control_name = row['control_name']
        iso_name = row['iso_name']
        pattern = row['pattern']
        results.append([ref_name, control_name, iso_name, pattern])

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
    pos_uniq_mut_aa = []
    for pos, mut_aa_stat in pos_mut_stat.items():
        mut_aa_stat = [(k, v) for k, v in mut_aa_stat.items()]
        if len(mut_aa_stat) == 1 and mut_aa_stat[0][-1] == len(results):
            pos_uniq_mut_aa.append(pos)
        most_common_mut_aa = sorted(mut_aa_stat, key=lambda x: x[1])[-1][0]
        pos_most_common_mut_aa[pos] = most_common_mut_aa

    headers = ['ref_name', 'control_name', 'iso_name']
    for pos in mut_pos_list:
        if pos in pos_uniq_mut_aa:
            continue

        ref_aa = pos_ref_aa_dict[pos]
        mut_aa = pos_most_common_mut_aa[pos]
        headers.append('{}{}{}'.format(ref_aa, pos, mut_aa))

    summary_results = []
    for ref_name, control_name, iso_name, pattern in results:
        result = {
            'ref_name': ref_name,
            'control_name': control_name,
            'iso_name': iso_name
        }

        mut_list = pattern.split('+')
        for mut in mut_list:
            ref_aa, pos, mut_aa = re.match(
                r'^(\D+)([\d\-]+)(\D+)$', mut).groups()
            if '-' in pos:
                pos = pos.split('-', 1)[0]

            pos = int(pos)
            if pos in pos_uniq_mut_aa:
                continue

            most_common_mut_aa = pos_most_common_mut_aa[pos]

            result['{}{}{}'.format(ref_aa, pos, most_common_mut_aa)] = (
                '-' if most_common_mut_aa == mut_aa else mut_aa
            )

        summary_results.append(result)

    summary_results.sort(key=lambda x: x['ref_name'])

    dump_csv(csv_save_path, summary_results, headers=headers)
