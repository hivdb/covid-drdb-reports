from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from preset import RESISTANCE_FILTER


MAIN_SQL = """
SELECT s.ref_name, s.rx_name, SUM(s.cumulative_count)
    FROM
    susc_results as s,
    {rxtype} AS rxtype
    WHERE rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
    {filters}
    GROUP BY s.ref_name, s.rx_name;
"""

ROWS = {
    '501Y': {
        'filter': [
            "AND s.strain_name = 'S:501Y'",
        ]
    },
    '∆69/70': {
        'filter': [
            "AND s.strain_name = 'S:69del+70del'",
        ]
    },
    '∆69/70 + 501Y': {
        'filter': [
            "AND s.strain_name = 'S:69del+70del+501Y'",
        ]
    },
    '∆69/70 + 501Y + 570D': {
        'filter': [
            "AND s.strain_name = 'S:69del+70del+501Y+570D'",
        ]
    },
    '∆69/70 + 501Y + 453F': {
        'filter': [
            "AND s.strain_name = 'S:69del+70del+453F'",
        ]
    },
    '484K': {
        'filter': [
            "AND s.strain_name = 'S:484K'"
        ]
    },
    '484K + 501Y': {
        'filter': [
            "AND s.strain_name = 'S:484K+501Y'"
        ]
    },
    '417N': {
        'filter': [
            "AND s.strain_name = 'S:417N'"
        ]
    },
    '417N + 484K + 501Y(B.1.351 RBD)': {
        'filter': [
            ("AND ("
             "   s.strain_name = 'S:417N+484K+501Y'"
             "   OR s.strain_name = 'B.1.351 RBD')"),
        ]
    },
    '439K': {
        'filter': [
            "AND s.strain_name = 'S:439K'"
        ]
    },
}

SUBROWS = {
    'CP': {
        'rxtype': 'rx_conv_plasma',
    },
    'IP': {
        'rxtype': 'rx_immu_plasma',
    },
}


def gen_tableS2(conn):
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
                    records.append({
                        'Strain name': row_name,
                        'CP name': i[1],
                        'Resistance level': resist_name,
                        '#Published': i[2],
                        'Reference': i[0]
                    })

    records.sort(key=itemgetter(
        'Strain name', 'CP name', 'Reference'))

    save_path = DATA_FILE_PATH / 'TableS2.csv'
    dump_csv(save_path, records)
