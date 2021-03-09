from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json
from preset import SYNONYM2AB_NAME
from preset import AB_NAME2MAB_CLASS
from operator import itemgetter
from collections import defaultdict
from preset import RESISTANCE_FILTER
from preset import MAB_RENAME
from preset import round_number
from preset import EXCLUDE_MAB


MAIN_SQL = """
SELECT  s.ref_name,
        s.rx_name,
        rxtype.class,
        s.fold_cmp,
        s.fold
    FROM
    susc_results as s,
    (
        SELECT DISTINCT _rxtype.ref_name, _rxtype.rx_name, ab.class
        FROM {rxtype} AS _rxtype,
        (
            SELECT a.ab_name, b.class, b.target, b.source FROM
                antibodies AS a LEFT JOIN
                antibody_targets AS b
                ON a.ab_name = b.ab_name
            WHERE a.availability IS NOT NULL
            OR a.pdb_id IS NOT NULL
        ) AS ab
        WHERE _rxtype.ab_name = ab.ab_name
    ) as rxtype
    WHERE rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
    AND s.control_strain_name IN ('Control', 'Wuhan', 'S:614G')
    -- AND s.ineffective IS NULL
    {filters}
    GROUP BY s.ref_name, s.rx_name;
"""

ROWS = {
    'N501Y': {
        'filter': [
            "AND s.strain_name = 'S:501Y'",
        ]
    },
    '∆69/70': {
        'filter': [
            "AND s.strain_name = 'S:69del+70del'",
        ]
    },
    '∆69/70 + N501Y': {
        'filter': [
            "AND s.strain_name = 'S:69del+70del+501Y'",
        ]
    },
    '∆69/70 + N501Y + A570D': {
        'filter': [
            "AND s.strain_name = 'S:69del+70del+501Y+570D'",
        ]
    },
    '∆69/70 + N501Y + Y453F': {
        'filter': [
            "AND s.strain_name = 'S:69del+70del+453F'",
        ]
    },
    '∆144': {
        'filter': [
            "AND s.strain_name = 'S:144del'",
        ]
    },
    'E484K': {
        'filter': [
            "AND s.strain_name = 'S:484K'"
        ]
    },
    'Y453F': {
        'filter': [
            "AND s.strain_name = 'S:453F'"
        ]
    },
    'L452R': {
        'filter': [
            "AND s.strain_name = 'S:452R'"
        ]
    },
    'E484K + N501Y': {
        'filter': [
            "AND s.strain_name = 'S:484K+501Y'"
        ]
    },
    'K417N': {
        'filter': [
            "AND s.strain_name = 'S:417N'"
        ]
    },
    'K417N + E484K + N501Y': {
        'filter': [
            ("AND ("
             "   s.strain_name = 'S:417N+484K+501Y'"
             "   OR s.strain_name = 'B.1.351 RBD')"),
        ]
    },
    'N439K': {
        'filter': [
            "AND s.strain_name = 'S:439K'"
        ]
    },
}


SUBROWS = {
    'mAb': {
        'rxtype': 'rx_antibodies',
    },
}


def gen_tableS7(conn):
    cursor = conn.cursor()

    records = []
    for row_name, attr_r in ROWS.items():
        for subrow_name, attr_subr in SUBROWS.items():
            for resist_name, resist_filter in RESISTANCE_FILTER.items():
                rxtype = attr_subr['rxtype']

                r_filter = attr_r.get('filter', [])
                filter = '\n    '.join(r_filter + resist_filter)
                sql = MAIN_SQL.format(
                    rxtype=rxtype,
                    filters=filter
                )
                # print(sql)

                cursor.execute(sql)
                for i in cursor.fetchall():
                    reference = i[0]

                    ab_name = SYNONYM2AB_NAME.get(i[1], i[1])
                    if ab_name in EXCLUDE_MAB:
                        continue

                    ab_class_info = AB_NAME2MAB_CLASS.get(ab_name)
                    ab_class = i[2]
                    if not ab_class and ab_class_info:
                        ab_class = ab_class_info['class']
                    if '/' in ab_name or '+' in ab_name:
                        ab_class = ''

                    fold = '{}'.format(round_number(i[4]))
                    records.append({
                        'Strain name': row_name,
                        'Mab name': MAB_RENAME.get(ab_name, ab_name),
                        'Class': ab_class or '',
                        # 'Resistance level': resist_name,
                        'Fold': fold,
                        'Reference': reference
                    })

    records.sort(key=itemgetter(
        'Strain name', 'Class', 'Mab name'))

    save_path = DATA_FILE_PATH / 'TableS7.csv'
    dump_csv(save_path, records)

    json_records = defaultdict(list)
    for r in records:
        strain = r['Strain name']
        json_records[strain].append({
            'strain': strain,
            'rx': r['Mab name'],
            'mab_class': r['Class'],
            'fold': r['Fold'].replace('>', '&gt;'),
            'reference': r['Reference']
        })

    records = []
    for strain, assays in json_records.items():
        records.append({
            'strain': strain,
            'assays': sorted(assays, key=itemgetter('mab_class')),
        })

    strain = sorted(records, key=itemgetter('strain'))
    save_path = DATA_FILE_PATH / 'tableS7.json'
    dump_json(save_path, records)
