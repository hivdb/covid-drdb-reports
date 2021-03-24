from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from collections import defaultdict
from preset import dump_json
from preset import RESISTANCE_FILTER
from preset import EXCLUDE_PLASMA
from preset import EXCLUDE_STUDIES
from preset import PLASMA_RENAME
from preset import PLASMA_POST_RENAME
from preset import RENAME_CP_EXECUTOR

from variant_filter import include_mutations


MAIN_SQL = """
SELECT s.ref_name, s.rx_name, SUM(s.cumulative_count), COUNT(1)
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
    'N501Y': {
        'filter': [
            include_mutations([
                'S:501Y',
                'S:501Y+614G'])
        ]
    },
    '∆69/70': {
        'filter': [
            include_mutations([
                'S:69del+70del',
                'S:69del+70del+614G'])
        ]
    },
    '∆69/70 + N501Y': {
        'filter': [
            include_mutations([
                'S:69del+70del+501Y',
                'S:69del+70del+501Y+614G'])
        ]
    },
    '∆69/70 + N501Y + A570D': {
        'filter': [
            include_mutations([
                'S:69del+70del+501Y+570D',
                'S:69del+70del+501Y+570D+614G'])
        ]
    },
    '∆69/70 + N501Y + Y453F': {
        'filter': [
            include_mutations([
                'S:69del+70del+453F',
                'S:69del+70del+453F+614G'])
        ]
    },
    '∆144': {
        'filter': [
            include_mutations([
                'S:144del',
                'S:144del+614G'])
        ]
    },
    'E484K': {
        'filter': [
            include_mutations([
                'S:484K',
                'S:484K+614G'])
        ]
    },
    'E484K + N501Y': {
        'filter': [
            include_mutations([
                'S:484K+501Y',
                'S:484K+501Y+614G'])
        ]
    },
    'Y453F': {
        'filter': [
            include_mutations([
                'S:453F',
                'S:453F+614G'])
        ]
    },
    'L452R': {
        'filter': [
            include_mutations([
                'S:452R',
                'S:452R+614G'])
        ]
    },
    'K417N': {
        'filter': [
            include_mutations([
                'S:417N',
                'S:417N+614G'])
        ]
    },
    'K417N + E484K + N501Y': {
        'filter': [
            include_mutations([
                'S:417N+484K+501Y',
                'S:417N+484K+501Y+614G',
                'B.1.351 RBD'
                ])
        ]
    },
    'N439K': {
        'filter': [
            include_mutations([
                'S:439K',
                'S:439K+614G'])
        ]
    },
}

SUBROWS = {
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
    'IP': {
        'rxtype': 'rx_immu_plasma',
    },
}


def gen_table_plasma_muts(conn):
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
                    num_results = i[2]
                    num_records = i[3]

                    aggregated_results = False
                    if num_records < num_results and num_results > 1:
                        aggregated_results = True

                    # if cp_name in EXCLUDE_PLASMA:
                    #     continue
                    # exclude_tester = EXCLUDE_STUDIES.get(reference)
                    # if exclude_tester and exclude_tester(cp_name):
                    #     continue

                    cp_name = PLASMA_RENAME.get(cp_name, cp_name)
                    rename_executors = RENAME_CP_EXECUTOR.get(reference, [])
                    for tester, new_name in rename_executors:
                        if tester(cp_name):
                            cp_name = new_name

                    key = '{}{}{}'.format(variant_name, cp_name, reference)

                    if aggregated_results:
                        reference = '{}†'.format(reference)

                    rec = records[key]
                    rec['Variant name'] = variant_name
                    rec['Plasma'] = PLASMA_POST_RENAME.get(cp_name, cp_name)
                    rec['S'] = rec.get('S', 0)
                    rec['I'] = rec.get('I', 0)
                    rec['R'] = rec.get('R', 0)
                    if resist_name == 'susceptible':
                        rec['S'] += num_results
                    elif resist_name == 'partial':
                        rec['I'] += num_results
                    else:
                        rec['R'] += num_results

                    rec['Reference'] = reference

    records = list(records.values())
    records.sort(key=itemgetter(
        'Variant name', 'Plasma', 'Reference'))

    save_path = DATA_FILE_PATH / 'table_plasma_muts.csv'
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
    save_path = DATA_FILE_PATH / 'table_plasma_muts.json'
    dump_json(save_path, records)
