from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import row2dict


SQL = """
SELECT
    DISTINCT
        susc.ref_name,
        susc.assay_name,
        art.year,
        art.doi,
        art.url
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
    AND
    susc.ref_name NOT LIKE 'Syed%'
;
"""


def gen_omicron_assay(
        conn,
        folder=DATA_FILE_PATH / 'omicron',
        file_name='omicron_assay_raw_data.csv'):
    cursor = conn.cursor()

    cursor.execute(SQL)

    table = row2dict(cursor.fetchall())

    dump_csv(folder / file_name, table)
