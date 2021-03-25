from preset import DATA_FILE_PATH
from operator import itemgetter
from preset import dump_csv
from preset import dump_json


ESCAPE_SCORE_SQL = """
SELECT *
FROM dms_ace2_binding
WHERE
gene = 'S'
{filter};
"""

ROWS = {
    'N501Y': {
        'filter': [
            "AND position = 501",
            "AND amino_acid = 'Y'",
        ]
    },
    'E484K': {
        'filter': [
            "AND position = 484",
            "AND amino_acid = 'K'",
        ]
    },
    'K417N': {
        'filter': [
            "AND position = 417",
            "AND amino_acid = 'N'",
        ]
    },
    'L452R': {
        'filter': [
            "AND position = 452",
            "AND amino_acid = 'R'",
        ]
    },
}


def gen_table_binding(
    conn,
    csv_save_path=DATA_FILE_PATH / 'table_binding.csv',
    json_save_path=DATA_FILE_PATH / 'table_binding.json'
        ):

    cursor = conn.cursor()

    results = []
    for mut, attr in ROWS.items():
        filter = attr.get('filter', [])
        filter = '\n    '.join(filter)
        sql = ESCAPE_SCORE_SQL.format(filter=filter)

        cursor.execute(sql)
        for row in cursor.fetchall():
            results.append({
                'mutation': mut,
                'ace2_binding': row['ace2_binding'],
                'expression': row['expression'],
                'ace2_contact': row['ace2_contact']
            })

    results.sort(key=itemgetter('mutation'))

    dump_csv(csv_save_path, results)

    dump_json(json_save_path, results)
