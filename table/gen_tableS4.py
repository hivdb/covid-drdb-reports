from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json
from operator import itemgetter
from preset import RESISTANCE_FILTER
from preset import EXCLUDE_PLASMA
from preset import PLASMA_RENAME
from preset import PLASMA_POST_RENAME
from collections import defaultdict
from preset import EXCLUDE_STUDIES
from preset import RENAME_CP_EXECUTOR


MAIN_SQL = """
SELECT s.ref_name, s.rx_name, SUM(s.cumulative_count)
    FROM
    susc_results as s,
    {rxtype} AS rxtype
    WHERE rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
    AND s.control_variant_name IN ('Control', 'Wuhan', 'S:614G')
    -- AND s.ineffective IS NULL
    {filters}
    GROUP BY s.ref_name, s.rx_name;
"""

ROWS = {
    'B.1.1.7': {
        'filter': [
            "AND s.variant_name = 'B.1.1.7 Spike'",
        ]
    },
    'B.1.1.7 authentic': {
        'filter': [
            "AND s.variant_name = 'B.1.1.7 authentic'",
        ]
    },
    'B.1.351': {
        'filter': [
            "AND s.variant_name = 'B.1.351 Spike'",
        ]
    },
    'B.1.351 authentic': {
        'filter': [
            "AND s.variant_name = 'B.1.351 authentic'",
        ]
    },
    'P.1': {
        'filter': [
            "AND s.variant_name = 'P.1 Spike'",
        ]
    },
    'P.1 authentic': {
        'filter': [
            "AND s.variant_name = 'P.1 authentic'",
        ]
    },
    'CAL.20C': {
        'filter': [
            "AND s.variant_name IN ("
            "    'B.1.427 authentic',"
            "    'B.1.429 authentic',"
            "    'B.1.429 Spike')",
        ]
    },
}

SUBROWS = {
    'CP': {
        'rxtype': 'rx_conv_plasma',
        'cp_filters': [
            "AND rxtype.infection IN ('Unknown', 'S:614G')",
        ]
    },
    'IP': {
        'rxtype': 'rx_immu_plasma',
    },
}


def gen_tableS4(conn):
    cursor = conn.cursor()

    records = defaultdict(dict)
    for row_name, attr_r in ROWS.items():
        for subrow_name, attr_subr in SUBROWS.items():
            for resist_name, resist_filter in RESISTANCE_FILTER.items():
                rxtype = attr_subr['rxtype']

                r_filter = attr_r.get('filter', [])
                filter = '\n    '.join(r_filter + resist_filter)

                if subrow_name.lower().startswith('cp'):
                    filter += '\n   '
                    filter += '\n   '.join(attr_subr.get('cp_filters', []))

                sql = MAIN_SQL.format(
                    rxtype=rxtype,
                    filters=filter
                )
                # print(sql)

                cursor.execute(sql)
                for i in cursor.fetchall():
                    variant_name = row_name
                    cp_name = i[1]
                    reference = i[0]

                    if cp_name in EXCLUDE_PLASMA:
                        continue
                    exclude_tester = EXCLUDE_STUDIES.get(reference)
                    if exclude_tester and exclude_tester(cp_name):
                        continue

                    cp_name = PLASMA_RENAME.get(cp_name, cp_name)
                    rename_executor = RENAME_CP_EXECUTOR.get(reference)
                    if rename_executor:
                        tester, new_name = rename_executor
                        if tester(cp_name):
                            cp_name = new_name

                    key = '{}{}{}'.format(variant_name, cp_name, reference)

                    if variant_name.endswith('authentic'):
                        reference = '{}*'.format(reference)
                        variant_name = variant_name.split()[0]

                    rec = records[key]
                    rec['Variant name'] = variant_name
                    rec['Plasma'] = PLASMA_POST_RENAME.get(cp_name, cp_name)
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
        'Variant name', 'Plasma', 'Reference'))

    save_path = DATA_FILE_PATH / 'TableS4.csv'
    dump_csv(save_path, records)

    json_records = defaultdict(list)
    for r in records:
        variant = r['Variant name']
        json_records[variant].append({
            'variant': variant,
            'rx': r['Plasma'],
            's_fold': r['S'],
            'i_fold': r['I'],
            'r_fold': r['R'],
            'reference': r['Reference']
        })

    records = []
    for variant, assays in json_records.items():
        records.append({
            'variant': variant,
            'assays': sorted(assays, key=itemgetter('rx')),
        })

    variant = sorted(records, key=itemgetter('variant'))
    save_path = DATA_FILE_PATH / 'tableS4.json'
    dump_json(save_path, records)
