from preset import DATA_FILE_PATH
from operator import itemgetter
from preset import dump_csv
from preset import dump_json
from collections import defaultdict


SQL = """
SELECT
    vaccine_name,
    var_name,
    design,
    efficacy,
    ref_name
FROM
    vaccine_efficacy
"""


def gen_vp_efficacy(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_vp_efficacy.csv'
        ):

    cursor = conn.cursor()
    cursor.execute(SQL)

    vaccine_groups = defaultdict(list)
    for rec in cursor.fetchall():
        vaccine = rec['vaccine_name']
        vaccine_groups[vaccine].append(rec)

    records = []
    for vaccine, rec_list in vaccine_groups.items():
        efficacy = [{
            'name': ('Pre-variant'
                     if rec['var_name'] == 'B.1' else rec['var_name']),
            'efficacy': '{} ({})'.format(
                rec['efficacy'],
                rec['design'],
            ),
            'ref_name': rec['ref_name'],
        } for rec in rec_list]

        records.append({
            'Vaccine': vaccine,
            'Variant': efficacy
        })

    dump_csv(csv_save_path, records)