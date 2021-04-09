from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json

from operator import itemgetter
from collections import defaultdict
from resistancy import RESISTANCE_FILTER
from preset import round_number

from .preset import SYNONYM2AB_NAME
from .preset import AB_NAME2MAB_CLASS
from .preset import MAB_RENAME
from .preset import EXCLUDE_MAB

from variant_filter import include_mutations


MAB_MUTS_SQL = """
SELECT  s.ref_name as ref_name,
        rxtype.ab_name as ab_name,
        rxtype.ab_class as ab_class,
        s.fold_cmp as fold_cmp,
        s.fold as fold,
        s.ineffective as ineffective
    FROM
    susc_results as s,
    (
        SELECT * FROM
        (
        SELECT
            _rxtype.ref_name,
            _rxtype.rx_name,
            group_concat(_rxtype.ab_name, "/") as ab_name,
            group_concat(ab.class, "+") as ab_class,
            group_concat(ab.pdb_id) as pdb_id,
            group_concat(ab.availability) as availability
        FROM {rxtype} AS _rxtype,
        (
            SELECT DISTINCT
                a.ab_name,
                a.pdb_id,
                a.availability,
                b.class,
                b.target
            FROM
                antibodies AS a LEFT JOIN
                antibody_targets AS b
                ON a.ab_name = b.ab_name
        ) AS ab
        WHERE _rxtype.ab_name = ab.ab_name
        GROUP BY _rxtype.ref_name, _rxtype.rx_name
        )
        WHERE
            availability IS NOT NULL
            OR pdb_id IS NOT NULL
    ) as rxtype
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
                'S:69del+70del+501Y+453F',
                'S:69del+70del+501Y+453F+614G'])
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
    'E484K + N501Y': {
        'filter': [
            include_mutations([
                'S:484K+501Y',
                'S:484K+501Y+614G'])
        ]
    },
    'K417N': {
        'filter': [
            include_mutations([
                'S:417N',
                'S:417N+614G'])
        ]
    },
    'F490S': {
        'filter': [
            include_mutations([
                'S:490S',
                'S:490S+614G'])
        ]
    },
    'S494P': {
        'filter': [
            include_mutations([
                'S:494P',
                'S:494P+614G'])
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
    'mAb': {
        'rxtype': 'rx_antibodies',
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
        for subrow_name, attr_subr in SUBROWS.items():
            for resist_name, resist_filter in RESISTANCE_FILTER.items():
                rxtype = attr_subr['rxtype']

                r_filter = attr_r.get('filter', [])
                filter = '\n    '.join(r_filter + resist_filter)
                sql = MAB_MUTS_SQL.format(
                    rxtype=rxtype,
                    filters=filter
                )
                # print(sql)

                cursor.execute(sql)
                for row in cursor.fetchall():
                    reference = row['ref_name']
                    print(row['ab_name'])
                    ab_name = SYNONYM2AB_NAME.get(row['ab_name'], row['ab_name'])
                    if ab_name in EXCLUDE_MAB:
                        continue

                    ab_class_info = AB_NAME2MAB_CLASS.get(ab_name)
                    ab_class = row['ab_class']
                    if not ab_class and ab_class_info:
                        ab_class = ab_class_info['class']
                    if '/' in ab_name or '+' in ab_name:
                        ab_class = ''

                    # ineffective = row['ineffective']
                    fold = row['fold']
                    # if ineffective:
                    #     fold = 100
                    fold = '{}'.format(round_number(fold))

                    records.append({
                        'Variant name': row_name,
                        'Mab name': MAB_RENAME.get(ab_name, ab_name),
                        'Class': ab_class or '',
                        # 'Resistance level': resist_name,
                        'Fold': fold,
                        'Reference': reference
                    })

    records.sort(key=itemgetter(
        'Variant name', 'Class', 'Mab name'))

    dump_csv(csv_save_path, records)

    json_records = defaultdict(list)
    for r in records:
        variant = r['Variant name']
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
