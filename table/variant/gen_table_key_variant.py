from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json
from collections import defaultdict

from variant_filter import include_similar_mutations
from variant_filter import exclude_similar_mutations
from variant.preset import CONTROL_VARIANTS_SQL
from variant.preset import get_iso_names_by_var_name


def get_iso_name_filter(var_name, selector='all'):

    iso_names = get_iso_names_by_var_name(var_name, selector)

    return "AND s.iso_name IN ({})".format(iso_names)


def get_iso_name_ignore(var_name, selector='all'):

    iso_names = get_iso_names_by_var_name(var_name, selector)

    return "AND s.iso_name NOT IN ({})".format(iso_names)


TABLE_SUMMARY_CP_SQL = """
SELECT
    SUM(s.cumulative_count)
FROM
    susc_results AS s,
    {rxtype} AS rxtype
    {joins}
    WHERE
    rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
    AND
    s.potency_type IN ('IC50', 'NT50')
    AND
    s.control_iso_name IN {control_variants}
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
    AND
    s.potency_type IN ('IC50', 'NT50')
    AND
    s.control_iso_name IN {control_variants}
    -- AND s.ineffective IS NULL
    {filters};
"""

TABLE_SUMMARY_ROWS = {
    'N501Y': {
        'filter': [
            include_similar_mutations([
                'S:501Y',
            ])
        ]
    },
    'E484K': {
        'filter': [
            include_similar_mutations([
                'S:484K',
            ])
        ]
    },
    'K417N': {
        'filter': [
            include_similar_mutations([
                'S:417N',
            ])
        ]
    },
    'K417T': {
        'filter': [
            include_similar_mutations([
                'S:417T',
            ])
        ]
    },
    'L452R': {
        'filter': [
            include_similar_mutations([
                'S:452R',
            ])
        ]
    },
    'Other individual mutations': {
        'join': [
            "isolates AS vs",
            (
                "(SELECT iso_name, COUNT(*) AS num_muts FROM "
                "isolate_mutations GROUP BY iso_name) AS sm"
            )
        ],
        'filter': [
            "AND vs.iso_name = s.iso_name",
            "AND vs.site_directed IS TRUE",
            "AND sm.num_muts = 1 AND sm.iso_name = s.iso_name",
            exclude_similar_mutations([
                'S:484K',
            ]),
            exclude_similar_mutations([
                'S:501Y',
            ]),
            exclude_similar_mutations([
                'S:417N',
            ]),
            exclude_similar_mutations([
                'S:417T',
            ]),
            exclude_similar_mutations([
                'S:452R',
            ]),
        ]
    },
    'All individual mutations': {
        'join': [
            "isolates AS vs",
            (
                "(SELECT iso_name, COUNT(*) AS num_muts FROM "
                "isolate_mutations GROUP BY iso_name) AS sm"
            )
        ],
        'filter': [
            "AND vs.iso_name = s.iso_name",
            "AND vs.site_directed IS TRUE",
            "AND sm.num_muts = 1 AND sm.iso_name = s.iso_name",
        ]
    },
    'B.1.1.7': {
        'filter': [
            get_iso_name_filter('Alpha')
        ]
    },
    'B.1.351': {
        'filter': [
            get_iso_name_filter('Beta')
        ]
    },
    'P.1': {
        'filter': [
            get_iso_name_filter('Gamma')
        ]
    },
    'B.1.427/9': {
        'filter': [
            get_iso_name_filter('Epsilon')
        ]
    },
    'B.1.526': {
        'filter': [
            get_iso_name_filter('Iota')
        ]
    },
    'Other combinations of mutations': {
        'join': [
            "isolates AS vs",
            (
                "(SELECT iso_name, COUNT(*) AS num_muts FROM "
                "isolate_mutations GROUP BY iso_name) AS sm"
            )
        ],
        'filter': [
            "AND vs.iso_name = s.iso_name",
            "AND sm.num_muts > 1 AND sm.iso_name = s.iso_name",
            get_iso_name_ignore('Alpha'),
            get_iso_name_ignore('Beta'),
            get_iso_name_ignore('Gamma'),
            get_iso_name_ignore('Epsilon'),
            get_iso_name_ignore('Iota'),
            exclude_similar_mutations([
                'S:484K',
            ]),
            exclude_similar_mutations([
                'S:501Y',
            ]),
        ]
    },
    "All combinations of mutations": {
        'join': [
            "isolates AS vs",
            (
                "(SELECT iso_name, COUNT(*) AS num_muts FROM "
                "isolate_mutations GROUP BY iso_name) AS sm"
            )
        ],
        'filter': [
            "AND vs.iso_name = s.iso_name",
            "AND sm.num_muts > 1 AND sm.iso_name = s.iso_name",
        ]
    }
}


TABLE_SUMMARY_COLUMNS = {
    'CP': {
        'rxtype': 'rx_conv_plasma',
        'cp_filters': [
            (
                "AND ("
                "      rxtype.infected_iso_name IN {}"
                "   OR rxtype.infected_iso_name IS NULL"
                "    )".format(CONTROL_VARIANTS_SQL)
            ),
        ]
    },
    'VP': {
        'rxtype': 'rx_vacc_plasma',
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
    }
}


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
                filters=filter,
                control_variants=CONTROL_VARIANTS_SQL,
            )
            if column_name.lower().startswith('mab'):
                abfilters = attr_c.get('ab_filters', [])
                abfilters = '\n     '.join(abfilters)

                sql = TABLE_SUMMARY_MAB_SQL.format(
                    rxtype=rxtype,
                    ab_filters=abfilters,
                    joins=join,
                    filters=filter,
                    control_variants=CONTROL_VARIANTS_SQL,
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
    save_path = DATA_FILE_PATH / 'summary_variant.csv'
    headers = [
        'Variant name',
        'CP',
        'VP',
        'mAbs phase3',
        'mAbs structure'
    ]
    dump_csv(save_path, variant_summary, headers)
