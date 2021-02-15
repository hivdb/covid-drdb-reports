from preset import DATA_FILE_PATH
from preset import dump_csv


TABLE3_MAIN_SQL = """
SELECT COUNT(*) FROM
    susc_results AS s,
    {rxtype} AS rxtype
    {joins}
    WHERE rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
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
            "AND s.strain_name = 'B.1.1.7 Spike'",
        ]
    },
    'B.1.351': {
        'filter': [
            "AND s.strain_name = 'B.1.351 Spike'",
        ]
    },
    'P.1': {
        'filter': [
            "AND s.strain_name = 'P.1 Spike'",
        ]
    },
    'B.1.1.7 + B.1.351': {
        'filter': [
            (
                "AND ("
                "      s.strain_name = 'B.1.1.7 Spike' "
                "   OR s.strain_name = 'B.1.351 Spike' "
                "   )"
            ),
        ]
    },
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
            "AND vs.site_directed IS TRUE",
            "AND sm.num_muts > 1 AND sm.strain_name = s.strain_name",
        ]
    }
}


TABLE3_COLUMNS = {
    'CP': {
        'rxtype': 'rx_conv_plasma',
    },
    'IP': {
        'rxtype': 'rx_immu_plasma',
    },
    'mAbs phase3': {
        'rxtype': 'rx_antibodies',
        'join': [
            "antibodies AS ab",
        ],
        'filter': [
            "AND rxtype.ab_name = ab.ab_name",
            "AND ab.availability = 'Phase 3'",
        ]
    },
    'mAbs structure': {
        'rxtype': 'rx_antibodies',
        'join': [
            "antibodies AS ab",
        ],
        'filter': [
            "AND rxtype.ab_name = ab.ab_name",
            "AND ab.pdb_id IS NOT NULL",
        ]
    }
}

FOOT_TABLE_ROWS = {
    '1 mutation': {
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
    'mutation combination': {
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
            "AND sm.num_muts > 1 AND sm.strain_name = s.strain_name",
        ]
    },
    'VOC': {
        'filter': [
            (
                "AND ("
                "      s.strain_name = 'B.1.1.7 Spike' "
                "   OR s.strain_name = 'B.1.351 Spike' "
                "   )"
            ),
        ]
    },
}

FOOT_TABLE_COLUMNS = {
    'mAbs without structure': {
        'rxtype': 'rx_antibodies',
        'join': [
            "antibodies AS ab",
        ],
        'filter': [
            "AND rxtype.ab_name = ab.ab_name",
            "AND ab.pdb_id IS NULL",
        ]
    },
    'CP': {
        'rxtype': 'rx_conv_plasma',
    },
    'VP': {
        'rxtype': 'rx_immu_plasma',
    },
}


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
            sql = TABLE3_MAIN_SQL.format(
                rxtype=rxtype,
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
                '#Published': result
            })

    for row_name, attr_r in FOOT_TABLE_ROWS.items():
        for column_name, attr_c in FOOT_TABLE_COLUMNS.items():
            rxtype = attr_c['rxtype']

            r_join = attr_r.get('join', [])
            c_join = attr_c.get('join', [])
            join = ',\n    '.join([''] + r_join + c_join)

            r_filter = attr_r.get('filter', [])
            c_filter = attr_c.get('filter', [])
            filter = '\n    '.join(r_filter + c_filter)
            sql = TABLE3_MAIN_SQL.format(
                rxtype=rxtype,
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
                '#Published': result
            })

    save_path = DATA_FILE_PATH / 'TableS3.csv'
    dump_csv(save_path, records)
