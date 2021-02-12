from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import SYNONYM2AB_NAME
from preset import AB_NAME2MAB_CLASS
from operator import itemgetter
from preset import RESISTANCE_FILTER
from preset import MAB_RENAME
from preset import round_number
from preset import EXCLUDE_MAB


MAIN_SQL = """
SELECT s.ref_name, s.rx_name, ab.class, ab.target, ab.source,
       s.fold_cmp, s.fold
    FROM
    susc_results as s,
    {rxtype} AS rxtype,
    (SELECT a.ab_name, b.class, b.target, b.source FROM
        antibodies AS a LEFT JOIN
        antibody_targets AS b
        ON a.ab_name = b.ab_name
        WHERE a.availability IS NOT NULL
        OR b.class IS NOT NULL
    ) AS ab
    WHERE rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
    AND rxtype.ab_name = ab.ab_name
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
    'E484K': {
        'filter': [
            "AND s.strain_name = 'S:484K'"
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
    'K417N + E484K + N501Y (B.1.351 RBD)': {
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


def gen_tableS6(conn):
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
                    ab_name = SYNONYM2AB_NAME.get(i[1], i[1])
                    if ab_name in EXCLUDE_MAB:
                        continue

                    ab_class_info = AB_NAME2MAB_CLASS.get(ab_name)
                    ab_class = i[2]
                    if not ab_class and ab_class_info:
                        ab_class = ab_class_info['class']
                    ab_target = i[3]
                    if not ab_target and ab_class_info:
                        ab_target = ab_class_info['target']
                    ab_source = i[4]
                    if not ab_source and ab_class_info:
                        ab_source = ab_class_info['source']
                    records.append({
                        'Strain name': row_name,
                        'Mab name': MAB_RENAME.get(ab_name, ab_name),
                        'Class': ab_class,
                        # 'Target': ab_target,
                        # 'Source': ab_source,
                        # 'Resistance level': resist_name,
                        'Fold': '{}'.format(round_number(i[6])),
                        'Reference': i[0]
                    })

    records.sort(key=itemgetter(
        'Strain name', 'Mab name', 'Reference'))

    save_path = DATA_FILE_PATH / 'TableS4.csv'
    dump_csv(save_path, records)
