from preset import group_records_by
from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import row2dict


SQL = """
SELECT
    DISTINCT
        susc.ref_name,
        art.year,
        art.doi,
        art.url,
        iso.var_name,
        mab.ab_name
FROM
    susc_results susc,
    articles art,
    isolates iso,
    rx_mab_view mab
WHERE
    susc.ref_name = art.ref_name
    AND
    susc.iso_name = iso.iso_name
    AND
    iso.var_name IN ('Omicron/BA.1', 'Omicron/BA.2', 'Omicron/BA.1.1')
    AND
    susc.ref_name = mab.ref_name
    AND
    susc.rx_name = mab.rx_name
    AND
    mab.availability IS NOT NULL
;
"""

AB_NAME_MAP = {
    'Bamlanivimab': 'BAM/ETE',
    'Etesevimab': 'BAM/ETE',
    'Bamlanivimab/Etesevimab': 'BAM/ETE',
    'Casirivimab': 'CAS/IMD',
    'Imdevimab': 'CAS/IMD',
    'Casirivimab/Imdevimab': 'CAS/IMD',
    'Cilgavimab': 'CIL/TIX',
    'Tixagevimab': 'CIL/TIX',
    'Cilgavimab/Tixagevimab': 'CIL/TIX',
    'Sotrovimab': 'SOT',
    'Regdanvimab': 'REG',
    'Adintrevimab': 'ADG20',
    'Bebtelovimab': 'BEB',
    'Amubarvimab': 'AMU/ROM',
    'Romlusevimab': 'AMU/ROM',
    'Amubarvimab/Romlusevimab': 'AMU/ROM',
    'C135': 'C135/C144',
    'C144': 'C135/C144',
    'C135/C144': 'C135/C144',
}


def gen_omicron_ref_info(
        conn,
        folder=DATA_FILE_PATH / 'omicron',
        file_name='omicron_ref_info_raw_data.csv'):
    cursor = conn.cursor()

    cursor.execute(SQL)

    table = row2dict(cursor.fetchall())

    dump_csv(folder / file_name, table)

    ba_1_list = set([
        r['ref_name']
        for r in table
        if r['var_name'] == 'Omicron/BA.1'
    ])
    ba_1_1_list = set([
        r['ref_name']
        for r in table
        if r['var_name'] == 'Omicron/BA.1.1'
    ])
    ba_2_list = set([
        r['ref_name']
        for r in table
        if r['var_name'] == 'Omicron/BA.2'
    ])

    ba_1_only_list = ba_1_list - (ba_2_list & ba_1_1_list)
    # ba_2_only_list = ba_2_list - (ba_1_list & ba_1_1_list)
    # ba_1_1_only_list = ba_2_list - (ba_2_list & ba_1_list)

    ba_1_and_ba_2_list = ba_1_list & ba_2_list
    ba_1_and_ba_1_1_list = ba_1_list & ba_1_1_list

    group_by_ref_name = group_records_by(table, 'ref_name')
    all_ab_name = sorted(list(set([
        r['ab_name']
        for r in table
    ])))

    detail_records = []
    for ref_name, ref_name_list in group_by_ref_name.items():
        rec = {
            'ref_name': ref_name,
            'doi': (
                ref_name_list[0]['doi']
                if ref_name_list[0]['doi'] else ref_name_list[0]['url']),
            'year': ref_name_list[0]['year'],
        }
        ref_ab_name = [
            r['ab_name']
            for r in ref_name_list
        ]

        rec['BA.1'] = 1 if ref_name in ba_1_list else 0
        rec['BA.2'] = 1 if ref_name in ba_2_list else 0
        rec['BA.1.1'] = 1 if ref_name in ba_1_1_list else 0
        for ab_name in all_ab_name:
            if ab_name not in AB_NAME_MAP.keys():
                continue
            if ab_name in ref_ab_name:
                rec[AB_NAME_MAP[ab_name]] = 1
            else:
                rec[AB_NAME_MAP[ab_name]] = 0

        detail_records.append(rec)

    dump_csv(folder / 'omicron_ref_info_detail.csv', detail_records)

    process_subvariant(
        folder,
        ba_1_list, ba_2_list, ba_1_1_list, ba_1_only_list,
        ba_1_and_ba_2_list, ba_1_and_ba_1_1_list)

    process_mab(
        table,
        exclude=lambda x: (x not in ba_1_list),
        file_name=folder / 'omicron_ref_info_BA_1_mab.csv')

    process_mab(
        table,
        exclude=lambda x: (x not in ba_2_list),
        file_name=folder / 'omicron_ref_info_BA_2_mab.csv')

    process_mab(
        table,
        exclude=lambda x: (x not in ba_1_1_list),
        file_name=folder / 'omicron_ref_info_BA_1_1_mab.csv')

    process_mab2(
        table,
        exclude=lambda x: (x not in ba_1_list),
        file_name=folder / 'omicron_ref_info_BA_1_mab2.csv')

    process_mab2(
        table,
        exclude=lambda x: (x not in ba_2_list),
        file_name=folder / 'omicron_ref_info_BA_2_mab2.csv')

    process_mab2(
        table,
        exclude=lambda x: (x not in ba_1_1_list),
        file_name=folder / 'omicron_ref_info_BA_1_1_mab2.csv')


def process_subvariant(
        folder,
        ba_1_list, ba_2_list, ba_1_1_list, ba_1_only_list,
        ba_1_and_ba_2_list, ba_1_and_ba_1_1_list):
    detail = []
    detail.append({
        'variant': 'BA.1',
        'num_ref': len(ba_1_list),
        'ref_names': ', '.join(sorted(list(ba_1_list))),
    })
    detail.append({
        'variant': 'BA.2',
        'num_ref': len(ba_2_list),
        'ref_names': ', '.join(sorted(list(ba_2_list))),
    })
    detail.append({
        'variant': 'BA.1.1',
        'num_ref': len(ba_1_1_list),
        'ref_names': ', '.join(sorted(list(ba_1_1_list))),
    })
    detail.append({
        'variant': 'BA.1 only',
        'num_ref': len(ba_1_only_list),
        'ref_names': ', '.join(sorted(list(ba_1_only_list))),
    })
    detail.append({
        'variant': 'BA.1 and BA.2',
        'num_ref': len(ba_1_and_ba_2_list),
        'ref_names': ', '.join(sorted(list(ba_1_and_ba_2_list))),
    })
    detail.append({
        'variant': 'BA.1 and BA.1.1',
        'num_ref': len(ba_1_and_ba_1_1_list),
        'ref_names': ', '.join(sorted(list(ba_1_and_ba_1_1_list))),
    })

    dump_csv(folder / 'omicron_ref_info_by_subvariant.csv', detail)


def process_mab(table, exclude=lambda x: False, file_name=None):

    detail_meta = []
    for ref_name, ref_name_list in group_records_by(table, 'ref_name').items():
        if exclude(ref_name):
            continue

        for rec in ref_name_list:
            ab_name = rec['ab_name']
            if ab_name not in AB_NAME_MAP:
                continue
            rec['_ab_name'] = AB_NAME_MAP[ab_name]
            detail_meta.append(rec)

    detail_ab_name = sorted(list(
        set([r['_ab_name'] for r in detail_meta])
    ))
    detail = []
    for ab_name in detail_ab_name:
        ref_list = set([
            r['ref_name'] for r in detail_meta
            if r['_ab_name'] == ab_name
            ])
        detail.append({
            '_ab_name': ab_name,
            'num_ref': len(ref_list),
            'ref_info': ', '.join(sorted(ref_list))
        })

    dump_csv(
        file_name,
        detail)


def process_mab2(table, exclude=lambda x: False, file_name=None):

    detail_meta = []
    for ref_name, ref_name_list in group_records_by(table, 'ref_name').items():
        if exclude(ref_name):
            continue

        for rec in ref_name_list:
            ab_name = rec['ab_name']
            if ab_name not in AB_NAME_MAP:
                continue
            rec['_ab_name'] = AB_NAME_MAP[ab_name]
            detail_meta.append(rec)

    detail_ab_name = sorted(list(
        set([r['_ab_name'] for r in detail_meta])
    ))

    detail = []
    for ref_name, ref_name_list in group_records_by(
            detail_meta, 'ref_name').items():

        new_rec = {
            'ref_name': ref_name,
        }
        for mab in detail_ab_name:
            ab_list = [
                i
                for i in ref_name_list
                if i['_ab_name'] == mab
            ]
            new_rec[mab] = 1 if len(ab_list) else 0
        detail.append(new_rec)

    dump_csv(
        file_name,
        detail)
