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

DOMAIN_PATTERN_LIST = {
    'RBD': [
        'E484K',
        'N501Y',
        'K417N',
        'L452R',
        'N439K',
        'S477N',
        'Y453F',
        'F490S',
        'E484Q',
    ],
    'NTD': [
        'HV69-70∆',
        'LAL242-244∆',
        'L18F',
        'D215G',
        'S247R',
        'R246I',
        'D80A',
        'Y144∆',
    ],
    'CTD': [
        'P681H',
        'A570D',
        'P681R',
        'H655Y',
        'Q675H',
    ],
    'S2': [
        'A701V',
        'S982A',
        'D1118H',
        'T716I',
        'D936Y',
    ]
}


def gen_table_variant_summary(conn):
    iso_type = 'isolate_mutations_single_s_mut_view'
    save_path = DATA_FILE_PATH / 'variant' / 'summary_single_cp.csv'
    by_single(conn, iso_type, save_path)

    iso_type = 'isolate_mutations_combo_s_mut_view'
    save_path = DATA_FILE_PATH / 'variant' / 'summary_combo_cp.csv'
    by_combo(conn, iso_type, save_path)


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
    domain_other_group = defaultdict(list)
    for mut_name, rx_list in mut_group.items():

        domain = rx_list[0]['domain']

        if mut_name not in DOMAIN_PATTERN_LIST[domain]:
            domain_other_group[domain].extend(rx_list)
            continue

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

    for domain, rx_list in domain_other_group.items():

        rx_group = defaultdict(int)
        for item in rx_list:
            rx = item['rx_name']
            num = item['num_fold']
            rx_group[rx] += num

        record = {
            'pattern': 'other',
            'position': 10000,
            'domain': domain,
            'num_ref_name': len(set(r['ref_name'] for r in rx_list)),
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

    save_path = DATA_FILE_PATH / 'variant' / 'summary_single.csv'
    dump_csv(save_path, save_results, headers)

    dump_json(DATA_FILE_PATH / 'variant' / 'summary_single.json', save_results)


def by_combo(conn, iso_type, save_path):

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
            mut_name = rec['pattern']
            rec = {i: rec[i] for i in rec.keys()}
            rec['rx_name'] = column_name
            mut_group[mut_name].append(rec)

    save_results = []
    for mut_name, record_list in mut_group.items():

        rx_group = defaultdict(int)
        for item in record_list:
            rx = item['rx_name']
            num = item['num_fold']
            rx_group[rx] += num

        record = {
            'pattern': record_list[0]['pattern'],
            'var_name': record_list[0]['var_name'] or '',
            'num_ref_name': len(
                set(r['ref_name'] for r in record_list)),
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
        'var_name',
        'cp',
        'vp',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        ))

    headers = [
        'pattern',
        'var_name',
        'num_ref_name',
        'cp',
        'vp',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        'all mAbs',
        'num_exp',
    ]
    save_path = DATA_FILE_PATH / 'variant' / 'summary_combo.csv'
    dump_csv(save_path, save_results, headers)

    var_name_group_combo = defaultdict(list)
    for record_list in mut_group.values():
        for rec in record_list:
            var_name = rec['var_name']

            if not var_name:
                continue

            var_name = var_name.split()[0]
            var_name = var_name.split('/')[0]
            var_name_group_combo[var_name].append(rec)

    merged_same_combo = []
    for var_name, record_list in var_name_group_combo.items():

        rx_group = defaultdict(int)
        for item in record_list:
            rx = item['rx_name']
            num = item['num_fold']
            rx_group[rx] += num

        rbd_muts = set()

        for rec in record_list:
            rbd_muts |= get_RBD_mutation(rec['pattern'])

        record = {
            'pattern': '+'.join(list(rbd_muts)),
            'var_name': var_name,
            'num_ref_name': len(set(r['ref_name'] for r in record_list)),
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
        merged_same_combo.append(record)

    merged_same_combo.sort(key=itemgetter(
        'var_name',
        'cp',
        'vp',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        ))

    save_path = DATA_FILE_PATH / 'variant' / 'summary_combo_by_var.csv'
    dump_csv(save_path, merged_same_combo, headers)


def get_RBD_mutation(mutation_list):
    rbd_muts = set()
    for mut in mutation_list.split('+'):
        pos_search = re.search(r'\d+', mut)
        if not pos_search:
            continue

        pos = int(pos_search.group())
        if pos in range(306, 535):
            rbd_muts.add(mut)

    return rbd_muts
