from preset import DATA_FILE_PATH
from preset import dump_csv
from collections import defaultdict
from operator import itemgetter
import re


TABLE_SUMMARY_SQL = """
SELECT
    s.iso_name,
    s.ref_name,
    s.cumulative_count num_fold,
    iso.*
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
    {filters}
;
"""


TABLE_SUMMARY_COLUMNS = {
    'cp': {
        'rx_type': 'rx_conv_plasma',
    },
    'vp': {
        'rx_type': 'rx_vacc_plasma',
    },
    'mAbs phase3': {
        'rx_type': 'rx_mab_view',
        'filters': [
            "AND rx.availability IS NOT NULL",
        ],
    },
    'mAbs structure': {
        'rx_type': 'rx_mab_view',
        'filters': [
            "AND rx.availability IS NULL",
            "AND rx.pdb_id IS NOT NULL"
        ],
    },
    'other mAbs': {
        'rx_type': 'rx_mab_view',
        'filters': [
            "AND rx.availability IS NULL",
            "AND rx.pdb_id IS NULL"
        ],
    },
}


def gen_table_single_mut_summary(conn):
    iso_type = 'isolate_mutations_single_s_mut_view'
    save_path = DATA_FILE_PATH / 'variant' / 'summary_single_cp.csv'
    by_single(conn, iso_type, save_path)


def by_single(conn, iso_type, save_path):

    cursor = conn.cursor()

    mut_group = defaultdict(list)

    for column_name, attr_c in TABLE_SUMMARY_COLUMNS.items():
        rx_type = attr_c['rx_type']

        c_filter = attr_c.get('filters', [])
        filter = '\n    '.join(c_filter)

        sql = TABLE_SUMMARY_SQL.format(
            rx_type=rx_type,
            filters=filter,
            iso_type=iso_type
        )

        cursor.execute(sql)
        for rec in cursor.fetchall():
            mut_name = rec['single_mut_name']
            rec = {i: rec[i] for i in rec.keys()}
            rec['rx_name'] = column_name
            rec['num_fold'] = rec['num_fold'] or 0
            mut_group[mut_name].append(rec)

    save_results = []
    for mut_name, rx_list in mut_group.items():

        rx_group = defaultdict(int)
        for item in rx_list:
            rx = item['rx_name']
            num = item['num_fold']
            rx_group[rx] += num

        record = {
            'pattern': mut_name,
            'ref': rx_list[0]['ref'],
            'position': rx_list[0]['position'],
            'aa': rx_list[0]['amino_acid'],
            'domain': rx_list[0]['domain'],
            'num_ref_name': len(set([
                r['ref_name']
                for r in rx_list
            ])),
            'cp': 0,
            'vp': 0,
            'mAbs phase3': 0,
            'mAbs structure': 0,
            'other mAbs': 0,
        }

        for rx, num in rx_group.items():
            record[rx] = num
        record['all mAbs'] = (
            record['mAbs phase3'] + record['mAbs structure'] +
            record['other mAbs'])
        record['num_exp'] = record['all mAbs'] + record['cp'] + record['vp']
        save_results.append(record)

    save_results.sort(key=itemgetter(
        'position',
        'cp',
        'vp',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        ))

    headers = [
        'pattern',
        'ref',
        'position',
        'aa',
        'domain',
        'num_ref_name',
        'cp',
        'vp',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        'all mAbs',
        'num_exp',
    ]

    save_path = DATA_FILE_PATH / 'variant' / 'summary_single_muts_all.csv'
    dump_csv(save_path, save_results, headers)
