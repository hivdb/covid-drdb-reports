from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json
from variant.preset import KEY_MUTATIONS

from operator import itemgetter
from collections import defaultdict
from resistancy import RESISTANCE_FILTER
from resistancy import round_fold


MAB_MUTS_SQL = """
SELECT
    s.ref_name,
    rx.ab_name,
    rx.class,
    s.fold_cmp,
    s.fold,
    s.ineffective
FROM
    susc_results_50_wt_view as s,
    rx_mab_view as rx,
    {iso_type} mut
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.iso_name = mut.iso_name
    AND
    s.fold IS NOT NULL
    AND (
        rx.availability IS NOT NULL
        OR rx.pdb_id IS NOT NULL
    )
    AND
    {filters}
"""


def gen_table_mab_muts(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_mab_muts.csv',
        ):
    cursor = conn.cursor()

    records = []
    for row_name, attr_r in KEY_MUTATIONS.items():
        for resist_name, resist_filter in RESISTANCE_FILTER.items():
            iso_type = attr_r.get('iso_type')

            r_filter = attr_r.get('filter', [])
            filter = '\n    '.join(r_filter + resist_filter)
            sql = MAB_MUTS_SQL.format(
                iso_type=iso_type,
                filters=filter,
            )
            # print(sql)

            cursor.execute(sql)
            for row in cursor.fetchall():
                ref_name = row['ref_name']
                ab_name = row['ab_name']
                ab_class = row['class']

                fold = row['fold']
                # ineffective = row['ineffective']
                # if ineffective:
                #     fold = 1000
                fold = '{}'.format(round_fold(fold))

                records.append({
                    'pattern': row_name,
                    'ab_name': ab_name,
                    'class': ab_class or '',
                    # 'Resistance level': resist_name,
                    'fold': fold,
                    'ref_name': ref_name
                })

    records.sort(key=itemgetter(
        'pattern',
        'class',
        'ab_name',
        ))

    dump_csv(csv_save_path, records)

    json_records = defaultdict(list)
    for r in records:
        variant = r['pattern']
        json_records[variant].append({
            'variant': variant,
            'rx_name': r['ab_name'],
            'mab_class': r['class'],
            # 'fold': r['fold'].replace('>', '&gt;'),
            'fold': r['fold'],
            'ref_name': r['ref_name']
        })

    records = []
    for pattern, assays in json_records.items():
        records.append({
            'pattern': pattern,
            'assays': sorted(assays, key=itemgetter('mab_class')),
        })

    records.sort(key=itemgetter('pattern'))
