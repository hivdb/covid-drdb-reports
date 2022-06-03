from preset import DATA_FILE_PATH
from preset import dump_csv
from sql import row2dict
from collections import defaultdict

SQL = """
SELECT
    *
FROM
    articles
WHERE
    ref_name IN (
        SELECT
            ref_name
        FROM
            susc_results
        UNION ALL
        SELECT
            ref_name
        FROM
            invivo_selection_results
        UNION ALL
        SELECT
            ref_name
        FROM
            invitro_selection_results
    )
;
"""


def is_preprint(doi):
    if doi.startswith('10.1101'):
        return True
    else:
        return False


def gen_ref_info(
        conn,
        save_path=DATA_FILE_PATH / 'reference' / 'ref_info.csv'):

    cursor = conn.cursor()

    cursor.execute(SQL)

    records = row2dict(cursor.fetchall())

    dump_csv(save_path, records)

    preprints = []
    for rec in records:
        doi = rec['doi']
        if not doi:
            continue
        if is_preprint(doi):
            preprints.append(rec)

    year_month_group = defaultdict(int)
    for rec in preprints:
        doi = rec['doi']
        date = doi.split('/', 1)[-1]
        year, month, day, _id = date.split('.')
        year_month_group['{}-{}'.format(year, month)] += 1

    summary = []
    for year_month, num_ref in year_month_group.items():
        summary.append({
            'year_month': year_month,
            'num_ref': num_ref
        })

    dump_csv(DATA_FILE_PATH / 'reference' / 'preprint_summary.csv', summary)
