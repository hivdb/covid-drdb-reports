from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json
from operator import itemgetter
from collections import defaultdict
from resistancy import RESISTANCE_FILTER
from resistancy import round_fold

from .preset import MAB_RENAME
from mab.preset import RX_MAB
from variant.preset import CONTROL_VARIANTS_SQL
from variant.preset import get_iso_names_by_var_name


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
    s.inhibition_pcnt != 90
    AND
    s.control_iso_name IN {control_variants}
    AND
    s.fold IS NOT NULL
    AND (
        rx.availability IS NOT NULL
        OR rx.pdb_id IS NOT NULL
    )
    {filters}
"""


def get_iso_name_filter(var_name, selector='all'):

    iso_names = get_iso_names_by_var_name(var_name, selector)

    return "AND s.iso_name in ({})".format(iso_names)


ROWS = {
    'B.1.1.7': {
        'filter': [
            get_iso_name_filter('Alpha', 'spike')
        ]
    },
    'B.1.1.7 full genome': {
        'filter': [
            get_iso_name_filter('Alpha', 'genome')
        ]
    },
    'B.1.351': {
        'filter': [
            get_iso_name_filter('Beta', 'spike')
        ]
    },
    'B.1.351 full genome': {
        'filter': [
            get_iso_name_filter('Beta', 'genome')
        ]
    },
    'P.1': {
        'filter': [
            get_iso_name_filter('Gamma', 'spike')
        ]
    },
    'P.1 full genome': {
        'filter': [
            get_iso_name_filter('Gamma', 'genome')
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
    'B.1.617.2': {
        'filter': [
            get_iso_name_filter('Delta')
        ]
    },
    'B.1.617.1': {
        'filter': [
            get_iso_name_filter('Kappa')
        ]
    },
}


def gen_table_mab_variant(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_mab_variant.csv',
        json_save_path=DATA_FILE_PATH / 'table_mab_variant.json'
        ):
    cursor = conn.cursor()

    records = []

    uniq_record_list = set()

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
                iso_name = row_name
                if iso_name.endswith('full genome'):
                    reference = '{}*'.format(reference)
                    iso_name = iso_name.split()[0]

                uniq_record = (
                    iso_name,
                    ab_name,
                    reference,
                    fold,
                )
                if uniq_record in uniq_record_list:
                    continue
                else:
                    uniq_record_list.add(uniq_record)

                records.append({
                    'Variant name': iso_name,
                    'Mab name': ab_name,
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
