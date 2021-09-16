from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json

from operator import itemgetter
from collections import defaultdict
from resistancy import RESISTANCE_FILTER
from resistancy import round_fold

from .preset import MAB_RENAME
from mab.preset import RX_MAB
from variant_filter import include_similar_mutations
from variant.preset import CONTROL_VARIANTS_SQL

MAB_MUTS_SQL = """
SELECT
    s.ref_name,
    rx.ab_name,
    rx.class,
    s.fold_cmp,
    s.fold,
    s.ineffective
FROM
    susc_results as s,
    ({rx_type}) as rx
ON
    s.ref_name = rx.ref_name AND
    s.rx_name = rx.rx_name
WHERE
    s.potency_type IN ('IC50', 'NT50')
    AND
    s.control_iso_name IN ({control_variants})
    AND s.fold IS NOT NULL
    {filters}
    AND (
        rx.availability IS NOT NULL
        OR rx.pdb_id IS NOT NULL
    )
"""

ROWS = {
    'N501Y': {
        'filter': [
            include_similar_mutations([
                'S:501Y',
            ])
        ]
    },
    '∆69/70': {
        'filter': [
            include_similar_mutations([
                'S:69del+70del',
            ])
        ]
    },
    '∆69/70 + N501Y': {
        'filter': [
            include_similar_mutations([
                'S:69del+70del+501Y',
            ])
        ]
    },
    '∆69/70 + N501Y + A570D': {
        'filter': [
            include_similar_mutations([
                'S:69del+70del+501Y+570D',
            ])
        ]
    },
    '∆69/70 + N501Y + Y453F': {
        'filter': [
            include_similar_mutations([
                'S:69del+70del+501Y+453F',
            ])
        ]
    },
    '∆144': {
        'filter': [
            include_similar_mutations([
                'S:144del',
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
    'Y453F': {
        'filter': [
            include_similar_mutations([
                'S:453F',
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
    'E484K + N501Y': {
        'filter': [
            include_similar_mutations([
                'S:484K+501Y',
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
    'F490S': {
        'filter': [
            include_similar_mutations([
                'S:490S',
            ])
        ]
    },
    'S494P': {
        'filter': [
            include_similar_mutations([
                'S:494P',
            ])
        ]
    },
    'K417N + E484K + N501Y': {
        'filter': [
            include_similar_mutations([
                'S:417N+484K+501Y',
                ])
        ]
    },
    'N439K': {
        'filter': [
            include_similar_mutations([
                'S:439K',
            ])
        ]
    },
    'T478K': {
        'filter': [
            include_similar_mutations([
                'S:478K',
            ])
        ]
    },
}


def gen_table_mab_muts(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_mab_muts.csv',
        json_save_path=DATA_FILE_PATH / 'table_mab_muts.json'
        ):
    cursor = conn.cursor()

    records = []
    for row_name, attr_r in ROWS.items():
        for resist_name, resist_filter in RESISTANCE_FILTER.items():
            r_filter = attr_r.get('filter', [])
            filter = '\n    '.join(r_filter + resist_filter)
            sql = MAB_MUTS_SQL.format(
                filters=filter,
                rx_type=RX_MAB,
                control_variants=CONTROL_VARIANTS_SQL,
            )
            # print(sql)

            cursor.execute(sql)
            for row in cursor.fetchall():
                reference = row['ref_name']
                ab_name = row['ab_name']
                ab_class = row['class']

                fold = row['fold']
                # ineffective = row['ineffective']
                # if ineffective:
                #     fold = 100
                fold = '{}'.format(round_fold(fold))

                ab_name = MAB_RENAME.get(ab_name, ab_name)
                records.append({
                    'pattern': row_name,
                    'Mab name': ab_name,
                    'Class': ab_class or '',
                    # 'Resistance level': resist_name,
                    'Fold': fold,
                    'Reference': reference
                })

    records.sort(key=itemgetter(
        'pattern', 'Class', 'Mab name'))

    dump_csv(csv_save_path, records)

    json_records = defaultdict(list)
    for r in records:
        variant = r['pattern']
        json_records[variant].append({
            'variant': variant,
            'rx': r['Mab name'],
            'mab_class': r['Class'],
            'fold': r['Fold'].replace('>', '&gt;'),
            'reference': r['Reference']
        })

    records = []
    for variant, assays in json_records.items():
        records.append({
            'variant': variant,
            'assays': sorted(assays, key=itemgetter('mab_class')),
        })

    variant = sorted(records, key=itemgetter('variant'))
    dump_json(json_save_path, records)
