from operator import itemgetter
from preset import DATA_FILE_PATH
from preset import dump_json
from sql import row2dict


SQL = """
    SELECT
        a.ab_name,
        GROUP_CONCAT(a.synonym, ',') as synonyms
    FROM
        antibody_synonyms a,
        antibodies b
    WHERE
        a.ab_name = b.ab_name
        AND
        b.availability IS NOT NULL
    GROUP BY
        a.ab_name
;
"""


def mab_synonyms(
        conn,
        json_save_path=DATA_FILE_PATH / 'mab' / 'mab_synonyms.json'):

    cursor = conn.cursor()
    cursor.execute(SQL)

    records = row2dict(cursor.fetchall())

    records.sort(key=itemgetter('ab_name'))

    dump_json(json_save_path, records)
