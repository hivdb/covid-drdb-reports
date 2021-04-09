from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json
from collections import defaultdict

from variant_filter import include_mutations
from variant_filter import exclude_mutations


TABLE_SUMMARY_CP_SQL = """
SELECT SUM(s.cumulative_count) FROM
-- SELECT COUNT(1) FROM
    susc_results AS s,
    {rxtype} AS rxtype
    {joins}
    WHERE rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
    AND s.control_variant_name IN ('Control', 'Wuhan', 'S:614G')
    -- AND s.ineffective IS NULL
    {filters};
"""

TABLE_SUMMARY_MAB_SQL = """
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
    AND s.control_variant_name IN ('Control', 'Wuhan', 'S:614G')
    -- AND s.ineffective IS NULL
    {filters};
"""

TABLE_SUMMARY_ROWS = {
    'N501Y': {
        'filter': [
            include_mutations([
                'S:501Y',
                'S:501Y+614G'])
        ]
    },
    'E484K': {
        'filter': [
            include_mutations([
                'S:484K',
                'S:484K+614G'])
        ]
    },
    'K417N': {
        'filter': [
            include_mutations([
                'S:417N',
                'S:417N+614G'])
        ]
    },
    'K417T': {
        'filter': [
            include_mutations([
                'S:417T',
                'S:417T+614G'])
        ]
    },
    'L452R': {
        'filter': [
            include_mutations([
                'S:452R',
                'S:452R+614G'])
        ]
    },
    'Other individual mutations': {
        'join': [
            "virus_variants AS vs",
            (
                "(SELECT variant_name, COUNT(*) AS num_muts FROM "
                "variant_mutations GROUP BY variant_name) AS sm"
            )
        ],
        'filter': [
            "AND vs.variant_name = s.variant_name",
            "AND vs.site_directed IS TRUE",
            "AND sm.num_muts = 1 AND sm.variant_name = s.variant_name",
            exclude_mutations([
                'S:484K',
                'S:484K+614G'
            ]),
            exclude_mutations([
                'S:501Y',
                'S:501Y+614G'
            ]),
            exclude_mutations([
                'S:417N',
                'S:417N+614G'
            ]),
            exclude_mutations([
                'S:417T',
                'S:417T+614G'
            ]),
            exclude_mutations([
                'S:452R',
                'S:452R+614G'
            ]),
        ]
    },
    'All individual mutations': {
        'join': [
            "virus_variants AS vs",
            (
                "(SELECT variant_name, COUNT(*) AS num_muts FROM "
                "variant_mutations GROUP BY variant_name) AS sm"
            )
        ],
        'filter': [
            "AND vs.variant_name = s.variant_name",
            "AND vs.site_directed IS TRUE",
            "AND sm.num_muts = 1 AND sm.variant_name = s.variant_name",
        ]
    },
    'B.1.1.7': {
        'filter': [
            "AND s.variant_name IN ('B.1.1.7 Spike', 'B.1.1.7 full genome')",
        ]
    },
    'B.1.351': {
        'filter': [
            "AND s.variant_name IN ('B.1.351 Spike', 'B.1.351 full genome')",
        ]
    },
    'P.1': {
        'filter': [
            "AND s.variant_name IN ('P.1 Spike', 'P.1 full genome')",
        ]
    },
    'B.1.427/9': {
        'filter': [
            "AND s.variant_name IN ("
            "    'B.1.427 full genome',"
            "    'B.1.429 full genome',"
            "    'B.1.429 Spike')",
        ]
    },
    'B.1.526': {
        'filter': [
            include_mutations([
                'B.1.526 Spike',
                'B.1.526 full genome',
            ])
        ]
    },
    # 'B.1.1.7 + B.1.351': {
    #     'filter': [
    #         (
    #             "AND ("
    #             "      s.variant_name IN ('B.1.1.7 Spike',
    #                                      'B.1.1.7 full genome')"
    #             "   OR s.variant_name IN ('B.1.351 Spike',
    #                                      'B.1.351 full genome')"
    #             "   OR s.variant_name IN ('P.1 Spike',
    #                                      'P.1 full genome')"
    #             "   )"
    #         ),
    #     ]
    # },
    'Other combinations of mutations': {
        'join': [
            "virus_variants AS vs",
            (
                "(SELECT variant_name, COUNT(*) AS num_muts FROM "
                "variant_mutations GROUP BY variant_name) AS sm"
            )
        ],
        'filter': [
            "AND vs.variant_name = s.variant_name",
            "AND sm.num_muts > 1 AND sm.variant_name = s.variant_name",
            exclude_mutations([
                'B.1.1.7 Spike',
                'B.1.1.7 full genome',
            ]),
            exclude_mutations([
                'B.1.351 Spike',
                'B.1.351 full genome',
            ]),
            exclude_mutations([
                'P.1 Spike',
                'P.1 full genome',
            ]),
            exclude_mutations([
                'B.1.427 full genome',
                'B.1.429 full genome',
                'B.1.429 Spike'
            ]),
            exclude_mutations([
                'S:484K',
                'S:484K+614G'
            ]),
            exclude_mutations([
                'S:501Y',
                'S:501Y+614G'
            ]),
            exclude_mutations([
                'B.1.526 Spike',
                'B.1.526 full genome',
            ])
        ]
    },
    "All combinations of mutations": {
        'join': [
            "virus_variants AS vs",
            (
                "(SELECT variant_name, COUNT(*) AS num_muts FROM "
                "variant_mutations GROUP BY variant_name) AS sm"
            )
        ],
        'filter': [
            "AND vs.variant_name = s.variant_name",
            "AND sm.num_muts > 1 AND sm.variant_name = s.variant_name",
        ]
    }
}


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
            "AND ab.pdb_id IS NOT NULL",
            "AND ab.availability IS NULL",
        ],
    }
}

# TABLE_SUMMARY_ADDITIONAL_ROWS = {
#     '1 mutation': {
#         'join': [
#             "virus_variants AS vs",
#             (
#                 "(SELECT variant_name, COUNT(*) AS num_muts FROM "
#                 "variant_mutations GROUP BY variant_name) AS sm"
#             )
#         ],
#         'filter': [
#             "AND vs.variant_name = s.variant_name",
#             "AND vs.site_directed IS TRUE",
#             "AND sm.num_muts = 1 AND sm.variant_name = s.variant_name",
#         ]
#     },
#     'mutation combination': {
#         'join': [
#             "virus_variants AS vs",
#             (
#                 "(SELECT variant_name, COUNT(*) AS num_muts FROM "
#                 "variant_mutations GROUP BY variant_name) AS sm"
#             )
#         ],
#         'filter': [
#             "AND vs.variant_name = s.variant_name",
#             "AND vs.site_directed IS TRUE",
#             "AND sm.num_muts > 1 AND sm.variant_name = s.variant_name",
#         ]
#     },
#     'VOC': {
#         'filter': [
#             (
#                 "AND ("
#                 "      s.variant_name IN ('B.1.1.7 Spike',
#                                          'B.1.1.7 full genome')"
#                 "   OR s.variant_name IN ('B.1.351 Spike'
#                                          'B.1.351 full genome')"
#                 "   OR s.variant_name IN ('P.1 Spike'
#                                          'P.1 full genome')"
#                 "   )"
#             ),
#         ]
#     },
# }

# TABLE_SUMMARY_ADDITIONAL_COLUMNS = {
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


def gen_table_key_variant(conn):
    cursor = conn.cursor()

    records = []

    for row_name, attr_r in TABLE_SUMMARY_ROWS.items():
        for column_name, attr_c in TABLE_SUMMARY_COLUMNS.items():
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
            result = cursor.fetchone()
            result = result[0]
            records.append({
                'Variant name': row_name,
                'Rx name': column_name,
                '#Published': result or 0
            })

    # for row_name, attr_r in TABLE_SUMMARY_ADDITIONAL_ROWS.items():
    #     for column_name, attr_c in TABLE_SUMMARY_ADDITIONAL_COLUMNS.items():
    #         rxtype = attr_c['rxtype']

    #         r_join = attr_r.get('join', [])
    #         c_join = attr_c.get('join', [])
    #         join = ',\n    '.join([''] + r_join + c_join)

    #         r_filter = attr_r.get('filter', [])
    #         c_filter = attr_c.get('filter', [])
    #         filter = '\n    '.join(r_filter + c_filter)

    #         sql = TABLE_SUMMARY_CP_SQL.format(
    #             rxtype=rxtype,
    #             joins=join,
    #             filters=filter
    #         )
    #         if column_name.lower().startswith('mab'):
    #             ab_filters = attr_c.get('ab_filters', [])
    #             ab_filters = '\n     '.join(ab_filters)

    #             sql = TABLE_SUMMARY_MAB_SQL.format(
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
    #             'Variant name': row_name,
    #             'Rx name': column_name,
    #             '#Published': result
    #         })

    save_path = DATA_FILE_PATH / 'table_summary.csv'
    dump_csv(save_path, records)

    json_info = defaultdict(dict)
    for item in records:
        variant = item['Variant name']
        rx = item['Rx name']
        if rx == 'mAbs phase3':
            rx = 'phase3'
        elif rx == 'mAbs structure':
            rx = 'structure'
        else:
            rx = rx.lower()
        count = item.get('#Published', 0)
        json_info[variant][rx] = count

    result = []
    for key, value in json_info.items():
        rec = {'variant': key}
        rec.update(value)
        result.append(rec)

    save_path = DATA_FILE_PATH / 'table_summary.json'
    dump_json(save_path, result)

    variant_summary = defaultdict(dict)
    for item in records:
        variant = item['Variant name']
        rx = item['Rx name']
        published = item['#Published']
        variant_summary[variant][rx] = published
        variant_summary[variant]['Variant name'] = variant

    variant_summary = list(variant_summary.values())
    save_path = DATA_FILE_PATH / 'table_summary_figure.csv'
    headers = [
        'Variant name',
        'CP',
        'VP',
        'mAbs phase3',
        'mAbs structure'
    ]
    dump_csv(save_path, variant_summary, headers)
