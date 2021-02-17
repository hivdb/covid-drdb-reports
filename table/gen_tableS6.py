from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from collections import defaultdict
from preset import dump_json
from preset import RESISTANCE_FILTER
from preset import EXCLUDE_PLASMA
from preset import PLASMA_RENAME
from preset import PLASMA_POST_RENAME


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
    'CP': {
        'rxtype': 'rx_conv_plasma',
    },
    'IP': {
        'rxtype': 'rx_immu_plasma',
    },
}


def gen_tableS6(conn):
    cursor = conn.cursor()

    records = defaultdict(dict)
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
                    strain_name = row_name
                    cp_name = i[1]
                    if cp_name in EXCLUDE_PLASMA:
                        continue
                    cp_name = PLASMA_RENAME.get(cp_name, cp_name)
                    reference = i[0]
                    key = '{}{}{}'.format(strain_name, cp_name, reference)

                    rec = records[key]
                    rec['Strain name'] = strain_name
                    rec['Plasma'] =  PLASMA_POST_RENAME.get(cp_name, cp_name)
                    rec['S'] = rec.get('S', 0)
                    rec['I'] = rec.get('I', 0)
                    rec['R'] = rec.get('R', 0)
                    if resist_name == 'susceptible':
                        rec['S'] += i[2]
                    elif resist_name == 'partial':
                        rec['I'] += i[2]
                    else:
                        rec['R'] += i[2]

                    rec['Reference'] = reference

    records = list(records.values())
    records.sort(key=itemgetter(
        'Strain name', 'Plasma', 'Reference'))

    save_path = DATA_FILE_PATH / 'TableS6.csv'
    dump_csv(save_path, records)

    json_records = defaultdict(list)
    for r in records:
        strain = r['Strain name']
        json_records[strain].append({
            'strain': strain,
            'rx': r['Plasma'],
            's_fold': r['S'],
            'i_fold': r['I'],
            'r_fold': r['R'],
            'reference': r['Reference']
        })

    records = []
    for strain, assays in json_records.items():
        records.append({
            'strain': strain,
            'assays': sorted(assays, key=itemgetter('rx')),
        })

    strain = sorted(records, key=itemgetter('strain'))
    save_path = DATA_FILE_PATH / 'tableS6.json'
    dump_json(save_path, records)
