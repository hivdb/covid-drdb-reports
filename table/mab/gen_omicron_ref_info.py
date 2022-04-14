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
    iso.var_name like 'Omicron%'
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
    'Amubarvimab': 'BRII-196/198',
    'Romlusevimab': 'BRII-196/198',
    'Amubarvimab/Romlusevimab': 'BRII-196/198',
    'C135': 'C135/C144',
    'C144': 'C135/C144',
    'C135/C144': 'C135/C144',
}


def gen_omicron_ref_info(
        conn,
        folder=DATA_FILE_PATH / 'mab',
        file_name='omicron_ref_info.csv'):
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

    ba_1_ref_info = set(ba_1_list) | set(ba_1_1_list)
    ba_1_only_ref_info = ba_1_ref_info - ba_2_list
    ba_2_only_ref_info = ba_2_list - ba_1_ref_info
    ba_1_and_ba_2_ref_info = ba_1_ref_info & ba_2_list

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

        rec['BA.1'] = 1 if ref_name in ba_1_ref_info else 0
        rec['BA.2'] = 1 if ref_name in ba_2_list else 0
        for ab_name in all_ab_name:
            if ab_name not in AB_NAME_MAP.keys():
                continue
            if ab_name in ref_ab_name:
                rec[AB_NAME_MAP[ab_name]] = 1
            else:
                rec[AB_NAME_MAP[ab_name]] = 0

        detail_records.append(rec)

    dump_csv(folder / 'omicron_ref_meta_info.csv', detail_records)

    detail = []
    detail.append({
        'variant': 'BA.1',
        'num_ref': len(ba_1_only_ref_info),
        'ref_names': ', '.join(sorted(list(ba_1_only_ref_info))),
    })
    detail.append({
        'variant': 'BA.2',
        'num_ref': len(ba_2_only_ref_info),
        'ref_names': ', '.join(sorted(list(ba_2_only_ref_info))),
    })
    detail.append({
        'variant': 'BA.1 and BA.2',
        'num_ref': len(ba_1_and_ba_2_ref_info),
        'ref_names': ', '.join(sorted(list(ba_1_and_ba_2_ref_info))),
    })

    dump_csv(folder / 'omicron_ref_info_about_variant.csv', detail)

    process_mab(
        table,
        filter=lambda x: (x not in ba_1_only_ref_info) or (
                x in ba_2_list),
        file_name=folder / 'omicron_ref_info_BA_1_mab.csv')

    process_mab(
        table,
        filter=lambda x: (x not in ba_2_list),
        file_name=folder / 'omicron_ref_info_BA_2_mab.csv')


def process_mab(table, filter=lambda x: False, file_name=None):

    detail_meta = []
    for ref_name, ref_name_list in group_records_by(table, 'ref_name').items():
        if filter(ref_name):
            continue

        for rec in ref_name_list:
            ab_name = rec['ab_name']
            if ab_name not in AB_NAME_MAP:
                continue
            rec['ab_name'] = AB_NAME_MAP[ab_name]
            detail_meta.append(rec)

    detail_ab_name = sorted(list(
        set([r['ab_name'] for r in detail_meta])
    ))
    detail = []
    for ab_name in detail_ab_name:
        ref_list = set([
            r['ref_name'] for r in detail_meta
            if r['ab_name'] == ab_name
            ])
        detail.append({
            'ab_name': ab_name,
            'num_ref': len(ref_list),
            'ref_info': ', '.join(sorted(ref_list))
        })

    dump_csv(
        file_name,
        detail)
