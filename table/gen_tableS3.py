from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json
from collections import defaultdict


TABLE3_MAIN_SQL = """
SELECT SUM(s.cumulative_count) FROM
-- SELECT COUNT(1) FROM
    susc_results AS s,
    {rxtype} AS rxtype
    {joins}
    WHERE rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
    AND s.control_strain_name IN ('Control', 'Wuhan', 'S:614G')
    -- AND s.ineffective IS NULL
    {filters};
"""

TABLE3_MAB_SQL = """
SELECT SUM(s.cumulative_count) FROM
-- SELECT COUNT(1) FROM
    susc_results AS s,
    (
        SELECT DISTINCT _rxtype.ref_name, _rxtype.rx_name
        FROM {rxtype} AS _rxtype, antibodies AS ab
        WHERE _rxtype.ab_name = ab.ab_name
        {ab_filters}
    ) as rxtype
    {joins}
    WHERE rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
    AND s.control_strain_name IN ('Control', 'Wuhan', 'S:614G')
    -- AND s.ineffective IS NULL
    {filters};
"""

TABLE3_ROWS = {
    'N501Y': {
        'filter': [
            "AND s.strain_name = 'S:501Y'"
        ]
    },
    'E484K': {
        'filter': [
            "AND s.strain_name = 'S:484K'"
        ]
    },
    'Other individual mutations': {
        'join': [
            "virus_strains AS vs",
            (
                "(SELECT strain_name, COUNT(*) AS num_muts FROM "
                "strain_mutations GROUP BY strain_name) AS sm"
            )
        ],
        'filter': [
            "AND vs.strain_name = s.strain_name",
            "AND vs.site_directed IS TRUE",
            "AND sm.num_muts = 1 AND sm.strain_name = s.strain_name",
            "AND s.strain_name != 'S:484K'",
            "AND s.strain_name != 'S:501Y'",
        ]
    },
    'All individual mutations': {
        'join': [
            "virus_strains AS vs",
            (
                "(SELECT strain_name, COUNT(*) AS num_muts FROM "
                "strain_mutations GROUP BY strain_name) AS sm"
            )
        ],
        'filter': [
            "AND vs.strain_name = s.strain_name",
            "AND vs.site_directed IS TRUE",
            "AND sm.num_muts = 1 AND sm.strain_name = s.strain_name",
        ]
    },
    'B.1.1.7': {
        'filter': [
            "AND s.strain_name IN ('B.1.1.7 Spike', 'B.1.1.7 authentic')",
        ]
    },
    'B.1.351': {
        'filter': [
            "AND s.strain_name IN ('B.1.351 Spike', 'B.1.351 authentic')",
        ]
    },
    'P.1': {
        'filter': [
            "AND s.strain_name IN ('P.1 Spike', 'P.1 authentic')",
        ]
    },
    'CAL.20C': {
        'filter': [
            "AND s.strain_name IN ("
            "    'B.1.427 authentic',"
            "    'B.1.429 authentic',"
            "    'B.1.429 Spike')",
        ]
    },
    # 'B.1.1.7 + B.1.351': {
    #     'filter': [
    #         (
    #             "AND ("
    #             "      s.strain_name IN ('B.1.1.7 Spike',
    #                                      'B.1.1.7 authentic')"
    #             "   OR s.strain_name IN ('B.1.351 Spike',
    #                                      'B.1.351 authentic')"
    #             "   OR s.strain_name IN ('P.1 Spike',
    #                                      'P.1 authentic')"
    #             "   )"
    #         ),
    #     ]
    # },
    'Other muation combinations': {
        'join': [
            "virus_strains AS vs",
            (
                "(SELECT strain_name, COUNT(*) AS num_muts FROM "
                "strain_mutations GROUP BY strain_name) AS sm"
            )
        ],
        'filter': [
            "AND vs.strain_name = s.strain_name",
            "AND sm.num_muts > 1 AND sm.strain_name = s.strain_name",
            "AND s.strain_name NOT IN ('B.1.1.7 Spike', 'B.1.1.7 authentic')",
            "AND s.strain_name NOT IN ('B.1.351 Spike', 'B.1.351 authentic')",
            "AND s.strain_name NOT IN ('P.1 Spike', 'P.1 authentic')",
            "AND s.strain_name NOT IN ("
            "    'B.1.427 authentic',"
            "    'B.1.429 authentic',"
            "    'B.1.429 Spike')",
        ]
    },
    "All combinations of mutations": {
        'join': [
            "virus_strains AS vs",
            (
                "(SELECT strain_name, COUNT(*) AS num_muts FROM "
                "strain_mutations GROUP BY strain_name) AS sm"
            )
        ],
        'filter': [
            "AND vs.strain_name = s.strain_name",
            "AND sm.num_muts > 1 AND sm.strain_name = s.strain_name",
        ]
    }
}


TABLE3_COLUMNS = {
    'CP': {
        'rxtype': 'rx_conv_plasma',
        'cp_filters': [
            "AND rxtype.variant IN ('Generic', 'S:614G')",
        ]
    },
    'IP': {
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
            "AND ab.pdb_id IS NOT NULL",
            "AND ab.availability IS NULL",
        ],
    }
}

# TABLES3_2_ROWS = {
#     '1 mutation': {
#         'join': [
#             "virus_strains AS vs",
#             (
#                 "(SELECT strain_name, COUNT(*) AS num_muts FROM "
#                 "strain_mutations GROUP BY strain_name) AS sm"
#             )
#         ],
#         'filter': [
#             "AND vs.strain_name = s.strain_name",
#             "AND vs.site_directed IS TRUE",
#             "AND sm.num_muts = 1 AND sm.strain_name = s.strain_name",
#         ]
#     },
#     'mutation combination': {
#         'join': [
#             "virus_strains AS vs",
#             (
#                 "(SELECT strain_name, COUNT(*) AS num_muts FROM "
#                 "strain_mutations GROUP BY strain_name) AS sm"
#             )
#         ],
#         'filter': [
#             "AND vs.strain_name = s.strain_name",
#             "AND vs.site_directed IS TRUE",
#             "AND sm.num_muts > 1 AND sm.strain_name = s.strain_name",
#         ]
#     },
#     'VOC': {
#         'filter': [
#             (
#                 "AND ("
#                 "      s.strain_name IN ('B.1.1.7 Spike',
#                                          'B.1.1.7 authentic')"
#                 "   OR s.strain_name IN ('B.1.351 Spike'
#                                          'B.1.351 authentic')"
#                 "   OR s.strain_name IN ('P.1 Spike'
#                                          'P.1 authentic')"
#                 "   )"
#             ),
#         ]
#     },
# }

# TABLES3_2_COLUMNS = {
#     'mAbs without structure': {
#         'rxtype': 'rx_antibodies',
#         'ab_filters': [
#             "AND ab.pdb_id IS NULL",
#         ],
#     },
#     'CP': {
#         'rxtype': 'rx_conv_plasma',
#     },
#     'VP': {
#         'rxtype': 'rx_immu_plasma',
#     },
# }


def gen_tableS3(conn):
    cursor = conn.cursor()

    records = []

    for row_name, attr_r in TABLE3_ROWS.items():
        for column_name, attr_c in TABLE3_COLUMNS.items():
            rxtype = attr_c['rxtype']

            r_join = attr_r.get('join', [])
            c_join = attr_c.get('join', [])
            join = ',\n    '.join([''] + r_join + c_join)

            r_filter = attr_r.get('filter', [])
            c_filter = attr_c.get('filter', [])
            filter = '\n    '.join(r_filter + c_filter)

            if column_name.lower().startswith('cp'):
                filter += '\n   '
                filter += '\n   '.join(attr_c.get('cp_filters', []))

            sql = TABLE3_MAIN_SQL.format(
                rxtype=rxtype,
                joins=join,
                filters=filter
            )
            if column_name.lower().startswith('mab'):
                abfilters = attr_c.get('ab_filters', [])
                abfilters = '\n     '.join(abfilters)

                sql = TABLE3_MAB_SQL.format(
                    rxtype=rxtype,
                    ab_filters=abfilters,
                    joins=join,
                    filters=filter
                )
            # print(sql)

            cursor.execute(sql)
            result = cursor.fetchone()
            result = result[0]
            records.append({
                'Strain name': row_name,
                'Rx name': column_name,
                '#Published': result or 0
            })

    # for row_name, attr_r in TABLES3_2_ROWS.items():
    #     for column_name, attr_c in TABLES3_2_COLUMNS.items():
    #         rxtype = attr_c['rxtype']

    #         r_join = attr_r.get('join', [])
    #         c_join = attr_c.get('join', [])
    #         join = ',\n    '.join([''] + r_join + c_join)

    #         r_filter = attr_r.get('filter', [])
    #         c_filter = attr_c.get('filter', [])
    #         filter = '\n    '.join(r_filter + c_filter)

    #         sql = TABLE3_MAIN_SQL.format(
    #             rxtype=rxtype,
    #             joins=join,
    #             filters=filter
    #         )
    #         if column_name.lower().startswith('mab'):
    #             ab_filters = attr_c.get('ab_filters', [])
    #             ab_filters = '\n     '.join(ab_filters)

    #             sql = TABLE3_MAB_SQL.format(
    #                 rxtype=rxtype,
    #                 ab_filters=ab_filters,
    #                 joins=join,
    #                 filters=filter,
    #             )

    #         # print(sql)

    #         cursor.execute(sql)
    #         result = cursor.fetchone()
    #         result = result[0]
    #         records.append({
    #             'Strain name': row_name,
    #             'Rx name': column_name,
    #             '#Published': result
    #         })

    save_path = DATA_FILE_PATH / 'TableS3.csv'
    dump_csv(save_path, records)

    json_info = defaultdict(dict)
    for item in records:
        strain = item['Strain name']
        rx = item['Rx name']
        if rx == 'mAbs phase3':
            rx = 'phase3'
        elif rx == 'mAbs structure':
            rx = 'structure'
        else:
            rx = rx.lower()
        count = item.get('#Published', 0)
        json_info[strain][rx] = count

    result = []
    for key, value in json_info.items():
        rec = {'strain': key}
        rec.update(value)
        result.append(rec)

    save_path = DATA_FILE_PATH / 'tableS3.json'
    dump_json(save_path, result)
