from preset import DATA_FILE_PATH
from preset import dump_csv
from collections import defaultdict
from operator import itemgetter

from .preset import INDIV_VARIANT
from .preset import MULTI_VARIANT


TABLE_SUMMARY_CP_SQL = """
SELECT
    variant_name,
    cumulative_count
FROM
    susc_results AS s,
    {rxtype} AS rxtype
    {joins}
    WHERE rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
    AND s.control_variant_name IN ('Control', 'Wuhan', 'S:614G')
    {filters};
"""

TABLE_SUMMARY_MAB_SQL = """
SELECT
    variant_name,
    cumulative_count
FROM
    susc_results AS s,
    (
        SELECT DISTINCT _rxtype.ref_name, _rxtype.rx_name
        FROM {rxtype} AS _rxtype, antibodies AS ab
        WHERE _rxtype.ab_name = ab.ab_name
        {ab_filters}
    ) as rxtype
    {joins}
    WHERE rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
    AND s.control_variant_name IN ('Control', 'Wuhan', 'S:614G')
    {filters};
"""


TABLE_SUMMARY_COLUMNS = {
    'CP': {
        'rxtype': 'rx_conv_plasma',
        'cp_filters': [
            (
                "AND ("
                "      rxtype.infection IN ('S:614G')"
                "   OR rxtype.infection IS NULL"
                "    )"
            ),
        ]
    },
    'VP': {
        'rxtype': 'rx_immu_plasma',
    },
    'mAbs phase3': {
        'rxtype': 'rx_antibodies',
        'ab_filters': [
            "AND ab.availability IS NOT NULL",
        ],
    },
    'mAbs structure': {
        'rxtype': 'rx_antibodies',
        'ab_filters': [
            (
                "AND ab.ab_name in " +
                "(SELECT ab_name FROM antibody_targets"
                " WHERE pdb_id IS NOT NULL)"),
            "AND ab.availability IS NULL",
        ],
    },
    'other mAbs': {
        'rxtype': 'rx_antibodies',
        'ab_filters': [
            (
                "AND ab.ab_name in " +
                "(SELECT ab_name FROM antibody_targets"
                " WHERE pdb_id IS NOT NULL)"),
            "AND ab.availability IS NULL",
        ],
    },
}


def gen_table_variant(conn):
    cursor = conn.cursor()

    indiv_results = defaultdict(list)
    multi_results = defaultdict(list)

    for column_name, attr_c in TABLE_SUMMARY_COLUMNS.items():
        rxtype = attr_c['rxtype']

        c_join = attr_c.get('join', [])
        join = ',\n    '.join([''] + c_join)

        c_filter = attr_c.get('filter', [])
        filter = '\n    '.join(c_filter)

        if column_name.lower().startswith('cp'):
            filter += '\n   '
            filter += '\n   '.join(attr_c.get('cp_filters', []))

        sql = TABLE_SUMMARY_CP_SQL.format(
            rxtype=rxtype,
            joins=join,
            filters=filter
        )
        if column_name.lower().startswith('mab'):
            abfilters = attr_c.get('ab_filters', [])
            abfilters = '\n     '.join(abfilters)

            sql = TABLE_SUMMARY_MAB_SQL.format(
                rxtype=rxtype,
                ab_filters=abfilters,
                joins=join,
                filters=filter
            )
        # print(sql)

        cursor.execute(sql)
        for rec in cursor.fetchall():
            variant = rec['variant_name']
            count_num = rec['cumulative_count']
            main_name = INDIV_VARIANT.get(variant)
            if main_name:
                disp_name = main_name['disp']
                indiv_results[disp_name].append({
                    'Variant name': disp_name,
                    'Rx name': column_name,
                    '#Published': count_num or 0
                })
            else:
                main_name = MULTI_VARIANT.get(variant)
                if not main_name:
                    continue
                disp_name = main_name['disp']
                multi_results[disp_name].append({
                    'Variant name': disp_name,
                    'Nickname': main_name['nickname'],
                    'Rx name': column_name,
                    '#Published': count_num or 0
                })

    # print(len(indiv_results))
    # print(len(multi_results))

    save_indiv = []
    for main_name, record_list in indiv_results.items():

        variant_info = INDIV_VARIANT[main_name]
        rx_group = defaultdict(int)
        for item in record_list:
            rx = item['Rx name']
            num = item['#Published']
            rx_group[rx] += num

        record = {
            'Variant name': main_name,
            'RefAA': variant_info['ref_aa'],
            'Position': variant_info['position'],
            'AA': variant_info['aa'],
            'Domain': variant_info['domain'],
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
        save_indiv.append(record)

    save_indiv.sort(key=itemgetter(
        'Position',
        'CP',
        'VP',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        ))

    save_multi = []
    for main_name, record_list in multi_results.items():
        rx_group = defaultdict(int)
        nickname = record_list[0]['Nickname']
        for item in record_list:
            rx = item['Rx name']
            num = item['#Published']
            rx_group[rx] += num

        record = {
            'Variant name': main_name,
            'Nickname': nickname,
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
        save_multi.append(record)

    save_multi.sort(key=itemgetter(
        'Nickname',
        'CP',
        'VP',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        ))

    headers = [
        'Variant name',
        'RefAA',
        'Position',
        'AA',
        'Domain',
        'CP',
        'VP',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        'all mAbs',
    ]

    save_path = DATA_FILE_PATH / 'table_variant_indiv_figure.csv'
    dump_csv(save_path, save_indiv, headers)

    headers = [
        'Variant name',
        'Nickname',
        'CP',
        'VP',
        'mAbs phase3',
        'mAbs structure',
        'other mAbs',
        'all mAbs',
    ]
    save_path = DATA_FILE_PATH / 'table_variant_multi_figure.csv'
    dump_csv(save_path, save_multi, headers)
