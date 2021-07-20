from preset import DATA_FILE_PATH
from operator import itemgetter
from preset import dump_csv
from preset import dump_json


SQL = """
SELECT
    vaccine_name,
    var_name,
    design,
    efficacy
FROM
    vaccine_efficacy
"""


def gen_vp_efficacy(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_vp_efficacy.csv',
        json_save_path=DATA_FILE_PATH / 'table_vp_efficacy.json',
        ):

    cursor = conn.cursor()
    cursor.execute(SQL)

    records = []
    for rec in cursor.fetchall():
        records.append({
            'Vaccine': rec['vaccine_name'],
            'Variant': rec['var_name'],
            'Efficacy': '{} ({})'.format(
                rec['efficacy'],
                rec['design'],
            )
        })

    records.sort(key=itemgetter(
        'Vaccine', 'Variant'))

    dump_csv(csv_save_path, records)
    dump_json(json_save_path, records)
