from collections import defaultdict
from preset import MAB_RENAME
from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json


MAB_SCORE_SQL = """
SELECT * FROM dms_escape_score;
"""

BINDING_SCORE_SQL = """
SELECT * FROM dms_ace2_binding;
"""


def gen_dms(conn,
            csv_save_path=DATA_FILE_PATH / 'table_dms.csv',
            json_save_path=DATA_FILE_PATH / 'table_dms.json'):
    cursor = conn.cursor()

    cursor.execute(BINDING_SCORE_SQL)

    dms_info = defaultdict(dict)
    for rec in cursor.fetchall():
        gene = rec['gene']
        position = rec['position']
        amino_acid = rec['amino_acid']
        mutation = '{}:{}{}'.format(gene, position, amino_acid)
        dms_info[mutation] = {
            'gene': gene,
            'position': position,
            'amino_acid': amino_acid,
            'ace2_binding': rec['ace2_binding'],
            'expression': rec['expression'],
            'ace2_contact': rec['ace2_contact'],
        }

    headers = set()

    cursor.execute(MAB_SCORE_SQL)
    for rec in cursor.fetchall():
        gene = rec['gene']
        position = rec['position']
        amino_acid = rec['amino_acid']
        mutation = '{}:{}{}'.format(gene, position, amino_acid)

        rx_name = rec['rx_name']
        rx_name = MAB_RENAME.get(rx_name, rx_name)
        headers.add(rx_name)

        escape_score = rec['escape_score']
        dms_info[mutation][rx_name] = escape_score

    headers = [
        'gene', 'position', 'amino_acid',
        'ace2_binding', 'expression', 'ace2_contact'] + sorted(list(headers))

    dump_csv(csv_save_path, list(dms_info.values()), headers)
    # dump_json(json_save_path, list(dms_info.values()))
