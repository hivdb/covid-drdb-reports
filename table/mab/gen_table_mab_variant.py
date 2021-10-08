from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json
from variant.preset import KEY_VARIANTS
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
    isolate_mutations_combo_s_mut_view mut
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





def gen_table_mab_variant(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_mab_variant.csv',
        json_save_path=DATA_FILE_PATH / 'table_mab_variant.json'
        ):
    cursor = conn.cursor()

    records = []

    uniq_record_list = set()

    for row_name, attr_r in KEY_VARIANTS.items():
        for resist_name, resist_filter in RESISTANCE_FILTER.items():
            r_filter = attr_r.get('filter', [])
            filter = '\n    '.join(r_filter + resist_filter)
            sql = MAB_MUTS_SQL.format(
                filters=filter,
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

                iso_name = row_name
                if 'full genome' in iso_name:
                    reference = '{}*'.format(reference)
                    iso_name = iso_name.split()[0]

                # uniq_record = (
                #     iso_name,
                #     ab_name,
                #     reference,
                #     fold,
                # )
                # if uniq_record in uniq_record_list:
                #     continue
                # else:
                #     uniq_record_list.add(uniq_record)

                records.append({
                    'pattern': iso_name,
                    'ab_name': ab_name,
                    'class': ab_class or '',
                    # 'Resistance level': resist_name,
                    'fold': fold,
                    'ref_name': reference
                })

    records.sort(key=itemgetter(
        'pattern',
        'class',
        'ab_name'))

    dump_csv(csv_save_path, records)

    json_records = defaultdict(list)
    for r in records:
        variant = r['pattern']
        json_records[variant].append({
            'variant': variant,
            'rx_name': r['ab_name'],
            'mab_class': r['class'],
            # 'fold': r['Fold'].replace('>', '&gt;'),
            'fold': r['fold'],
            'ref_name': r['ref_name']
        })

    records = []
    for variant, assays in json_records.items():
        records.append({
            'variant': variant,
            'assays': sorted(assays, key=itemgetter('mab_class')),
        })

    variant = sorted(records, key=itemgetter('variant'))
    dump_json(json_save_path, records)
