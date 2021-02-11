from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from preset import RESISTANCE_FILTER
from preset import EXCLUDE_PLASMA
from preset import PLASMA_RENAME
from collections import defaultdict


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
    'B.1.1.7': {
        'filter': [
            ("AND (s.strain_name = 'B.1.1.7 Spike'"
                #  "  OR s.strain_name = 'B.1.1.7 S1'"
             ")"),
        ]
    },
    'B.1.351': {
        'filter': [
            "AND s.strain_name = 'B.1.351 Spike'",
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


def gen_tableS1(conn):
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
                    rec['Plasma'] = cp_name
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

    save_path = DATA_FILE_PATH / 'TableS1.csv'
    dump_csv(save_path, records)