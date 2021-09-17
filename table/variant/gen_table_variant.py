from preset import DATA_FILE_PATH
from preset import dump_csv
from collections import defaultdict
from operator import itemgetter

from .preset import ONE_MUT_VARIANT
from .preset import COMBO_MUT_VARIANT
from mab.preset import RX_MAB
from variant.preset import CONTROL_VARIANTS_SQL


TABLE_SUMMARY_CP_SQL = """
SELECT
    s.iso_name,
    s.ref_name,
    s.cumulative_count
FROM
    susc_results AS s,
    {rxtype} AS rxtype
ON
    rxtype.ref_name = s.ref_name
    AND rxtype.rx_name = s.rx_name
;
"""
# WHERE
#     s.potency_type IN ('IC50', 'NT50')
#     AND
#     s.control_iso_name IN ({control_variants})
#     AND
#     s.fold IS NOT NULL
#     {filters}


TABLE_SUMMARY_MAB_SQL = """
SELECT
    s.iso_name,
    s.ref_name,
    s.cumulative_count
FROM
    susc_results AS s,
    ({rxtype}) as rx
ON
    s.ref_name = rx.ref_name
    AND s.rx_name = rx.rx_name
WHERE
    1
    {filters}
;
"""


TABLE_SUMMARY_COLUMNS = {
    'CP': {
        'rxtype': 'rx_conv_plasma',
    },
    'VP': {
        'rxtype': 'rx_vacc_plasma',
    },
    'mAbs phase3': {
        'filters': [
            "AND rx.availability IS NOT NULL",
        ],
    },
    'mAbs structure': {
        'filters': [
            "AND rx.availability IS NULL",
            "AND rx.pdb_id IS NOT NULL"
        ],
    },
    'other mAbs': {
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
        '69-70∆',
        '242-244∆',
        'L18F',
        'D215G',
        'S247R',
        'R246I',
        'D80A',
        '141-145∆',
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


def gen_table_variant(conn):
    cursor = conn.cursor()

    indiv_results = defaultdict(list)
    combo_results = defaultdict(list)

    for column_name, attr_c in TABLE_SUMMARY_COLUMNS.items():
        c_join = attr_c.get('join', [])
        join = ',\n    '.join([''] + c_join)

        c_filter = attr_c.get('filters', [])
        filter = '\n    '.join(c_filter)

        if column_name.lower() in ['cp', 'vp']:
            rxtype = attr_c['rxtype']
            filter += '\n   '
            filter += '\n   '.join(attr_c.get('cp_filters', []))

            sql = TABLE_SUMMARY_CP_SQL.format(
                rxtype=rxtype,
                joins=join,
                filters=filter,
                control_variants=CONTROL_VARIANTS_SQL,
            )
        if 'mab' in column_name.lower():
            sql = TABLE_SUMMARY_MAB_SQL.format(
                rxtype=RX_MAB,
                filters=filter,
                control_variants=CONTROL_VARIANTS_SQL,
            )
        # print(sql)

        cursor.execute(sql)
        for rec in cursor.fetchall():
            variant = rec['iso_name']
            ref_name = rec['ref_name']
            num_results = rec['cumulative_count']
            main_name = ONE_MUT_VARIANT.get(variant)
            if main_name:
                disp_name = main_name['disp']
                indiv_results[disp_name].append({
                    'pattern': disp_name,
                    'rx_name': column_name,
                    'ref_name': ref_name,
                    'num_results': num_results or 0
                })
            else:
                main_name = COMBO_MUT_VARIANT.get(variant)
                if not main_name:
                    continue
                disp_name = main_name['disp']
                combo_results[disp_name].append({
                    'pattern': disp_name,
                    'var_name': main_name['var_name'],
                    'rx_name': column_name,
                    'ref_name': ref_name,
                    'num_results': num_results or 0
                })

    # print(len(indiv_results))
    # print(len(combo_results))

    save_indiv = []
    domain_other_group = defaultdict(list)
    for main_name, record_list in indiv_results.items():
        variant_info = ONE_MUT_VARIANT[main_name]

        domain = variant_info['domain']
        if main_name not in DOMAIN_PATTERN_LIST[domain]:
            domain_other_group[domain].extend(record_list)
            continue

        rx_group = defaultdict(int)
        for item in record_list:
            rx = item['rx_name']
            num = item['num_results']
            rx_group[rx] += num

        record = {
            'pattern': main_name,
            'RefAA': variant_info['ref_aa'],
            'Position': variant_info['position'],
            'AA': variant_info['aa'],
            'Domain': domain,
            'num_ref_name': len(set(r['ref_name'] for r in record_list)),
            'CP': 0,
            'VP': 0,
            'mAbs phase3': 0,
            'mAbs structure': 0,
            'other mAbs': 0,
        }
        for rx, num in rx_group.items():
            record[rx] = num
        record['all mAbs'] = (
            record['mAbs phase3'] + record['mAbs structure'] +
            record['other mAbs'])
        record['num_exp'] = record['all mAbs'] + record['CP'] + record['VP']
        save_indiv.append(record)

    for domain, record_list in domain_other_group.items():

        rx_group = defaultdict(int)
        for item in record_list:
            rx = item['rx_name']
            num = item['num_results']
            rx_group[rx] += num

        record = {
            'pattern': 'other',
            'Position': 10000,
            'Domain': domain,
            'num_ref_name': len(set(r['ref_name'] for r in record_list)),
            'CP': 0,
            'VP': 0,
            'mAbs phase3': 0,
            'mAbs structure': 0,
            'other mAbs': 0,
        }
        for rx, num in rx_group.items():
            record[rx] = num
        record['all mAbs'] = (
            record['mAbs phase3'] + record['mAbs structure'] +
            record['other mAbs'])
        record['num_exp'] = record['all mAbs'] + record['CP'] + record['VP']
        save_indiv.append(record)

    save_indiv.sort(key=itemgetter(
        'Position',
        'CP',
        'VP',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        ))

    headers = [
        'pattern',
        'RefAA',
        'Position',
        'AA',
        'Domain',
        'num_ref_name',
        'CP',
        'VP',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        'all mAbs',
        'num_exp',
    ]

    save_path = DATA_FILE_PATH / 'summary_variant_indiv.csv'
    dump_csv(save_path, save_indiv, headers)

    save_combo = []
    for main_name, record_list in combo_results.items():
        rx_group = defaultdict(int)
        var_name = record_list[0]['var_name']
        for item in record_list:
            rx = item['rx_name']
            num = item['num_results']
            rx_group[rx] += num

        record = {
            'pattern': main_name,
            'var_name': var_name,
            'num_ref_name': len(set(r['ref_name'] for r in record_list)),
            'CP': 0,
            'VP': 0,
            'mAbs phase3': 0,
            'mAbs structure': 0,
            'other mAbs': 0,
        }
        for rx, num in rx_group.items():
            record[rx] = num
        record['all mAbs'] = (
            record['mAbs phase3'] + record['mAbs structure'] +
            record['other mAbs'])
        record['num_exp'] = record['all mAbs'] + record['CP'] + record['VP']
        save_combo.append(record)

    save_combo.sort(key=itemgetter(
        'var_name',
        'CP',
        'VP',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        ))

    headers = [
        'pattern',
        'var_name',
        'num_ref_name',
        'CP',
        'VP',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        'all mAbs',
        'num_exp',
    ]
    save_path = DATA_FILE_PATH / 'summary_variant_combo.csv'
    dump_csv(save_path, save_combo, headers)

    var_name_group_combo = defaultdict(list)
    for record_list in combo_results.values():
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
            num = item['num_results']
            rx_group[rx] += num

        record = {
            'var_name': var_name,
            'num_ref_name': len(set(r['ref_name'] for r in record_list)),
            'CP': 0,
            'VP': 0,
            'mAbs phase3': 0,
            'mAbs structure': 0,
            'other mAbs': 0,
        }
        for rx, num in rx_group.items():
            record[rx] = num
        record['all mAbs'] = (
            record['mAbs phase3'] + record['mAbs structure'] +
            record['other mAbs'])
        record['num_exp'] = record['all mAbs'] + record['CP'] + record['VP']
        merged_same_combo.append(record)

    merged_same_combo.sort(key=itemgetter(
        'var_name',
        'CP',
        'VP',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        ))

    save_path = DATA_FILE_PATH / 'summary_variant_combo_by_var.csv'
    dump_csv(save_path, merged_same_combo, headers)
