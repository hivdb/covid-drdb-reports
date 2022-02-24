from preset import group_records_by
from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import row2dict


SQL = """
SELECT
    DISTINCT
        susc.ref_name,
        iso.var_name,
        mab.ab_name
FROM
    susc_results susc,
    isolates iso,
    rx_mab_view mab
WHERE
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


def gen_omicron_ref_info(
        conn,
        folder=DATA_FILE_PATH / 'mab',
        file_name='omicron_ref_info.csv'):
    cursor = conn.cursor()

    cursor.execute(SQL)

    table = row2dict(cursor.fetchall())

    dump_csv(folder / file_name, table)

    ba_1_list = [
        r['ref_name']
        for r in table
        if r['var_name'] == 'Omicron/BA.1'
    ]
    ba_1_1_list = [
        r['ref_name']
        for r in table
        if r['var_name'] == 'Omicron/BA.1+R346K'
    ]
    ba_2_list = [
        r['ref_name']
        for r in table
        if r['var_name'] == 'Omicron/BA.2'
    ]
    ref_names = sorted(list(
        set(ba_1_list) | set(ba_1_list) | set(ba_2_list))
    )

    detail_1 = []
    for ref_name in ref_names:
        rec = {
            'ref_name': ref_name,
            'BA.1': True if ref_name in ba_1_list else False,
            'BA.1.1': True if ref_name in ba_1_1_list else False,
            'BA.2': True if ref_name in ba_2_list else False,
        }
        detail_1.append(rec)

    dump_csv(folder / 'omicron_ref_info_detail_1.csv', detail_1)

    ab_list = sorted(list(set([
        r['ab_name']
        for r in table
    ])))

    detail_2 = []
    for ref_name, ref_name_list in group_records_by(table, 'ref_name').items():
        if (ref_name not in ba_1_1_list) and (ref_name not in ba_1_list) or (
            ref_name in ba_2_list):
            continue

        ab_name_list = [
            r['ab_name']
            for r in ref_name_list
        ]
        rec = {
            'ref_name': ref_name
        }
        for ab_name in ab_list:
            if ab_name in ab_name_list:
                rec[ab_name] = True
            else:
                rec[ab_name] = False
        detail_2.append(rec)

    dump_csv(
        folder / 'omicron_ref_info_detail_2.csv', detail_2,
        headers=['ref_name'] + ab_list)

    detail_3 = []
    for ref_name, ref_name_list in group_records_by(table, 'ref_name').items():
        if (ref_name not in ba_2_list):
            continue

        ab_name_list = [
            r['ab_name']
            for r in ref_name_list
        ]
        rec = {
            'ref_name': ref_name
        }
        for ab_name in ab_list:
            if ab_name in ab_name_list:
                rec[ab_name] = True
            else:
                rec[ab_name] = False
        detail_3.append(rec)

    dump_csv(
        folder / 'omicron_ref_info_detail_3.csv', detail_3,
        headers=['ref_name'] + ab_list)
