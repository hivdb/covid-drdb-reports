from preset import DATA_FILE_PATH
from operator import itemgetter
from preset import dump_csv
from preset import dump_json

from .preset import SYNONYM2AB_NAME
from .preset import MAB_RENAME


EUA_RX_SQL = """
SELECT
    a.rx_name
FROM
    rx_dms AS a,
    antibodies AS b
ON a.ab_name = b.ab_name
WHERE b.availability = 'Phase 3'
"""

ESCAPE_SCORE_SQL = """
SELECT *
FROM dms_escape_results
WHERE rx_name IN (""" + EUA_RX_SQL + """)
    AND gene = 'S'
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


def gen_escape(
    conn,
    csv_save_path=DATA_FILE_PATH / 'table_escape.csv',
    json_save_path=DATA_FILE_PATH / 'table_escape.json'
        ):

    cursor = conn.cursor()

    results = []
    for mut, attr in ROWS.items():
        filter = attr.get('filter', [])
        filter = '\n    '.join(filter)
        sql = ESCAPE_SCORE_SQL.format(filter=filter)

        cursor.execute(sql)
        for row in cursor.fetchall():
            rx_name = row['rx_name']
            rx_name = SYNONYM2AB_NAME.get(rx_name, rx_name)
            rx_name = MAB_RENAME.get(rx_name, rx_name)
            results.append({
                'mutation': mut,
                'rx_name': rx_name,
                'escape_score': row['escape_score']
            })

    results.sort(key=itemgetter('mutation', 'rx_name'))

    dump_csv(csv_save_path, results)

    dump_json(json_save_path, results)
