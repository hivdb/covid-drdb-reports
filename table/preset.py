import csv
from pathlib import Path

WS = Path(__file__).absolute().parent.parent
DATA_FILE_PATH = WS / 'susceptibility-data_files' / 'table'
DATA_FILE_PATH.mkdir(exist_ok=True)


def dump_csv(file_path, records, headers=[]):
    if not records:
        return
    if not headers and records:
        headers = records[0].keys()

    with open(file_path, 'w', encoding='utf-8-sig') as fd:
        writer = csv.DictWriter(fd, fieldnames=headers)
        writer.writeheader()
        writer.writerows(records)


SYNONYM2AB_NAME = {}
SYNONYM_GROUP = {}


def init_synonyms_map(conn):
    cursor = conn.cursor()
    sql = """
        SELECT * FROM antibody_synonyms
    """
    cursor.execute(sql)

    for item in cursor.fetchall():
        ab_name = item[0]
        synonym = item[1]
        SYNONYM2AB_NAME[synonym] = ab_name

        if synonym in SYNONYM_GROUP:
            SYNONYM_GROUP[synonym].add(ab_name)
            SYNONYM_GROUP[ab_name] = SYNONYM_GROUP[synonym]
        elif ab_name in SYNONYM_GROUP:
            SYNONYM_GROUP[ab_name].add(synonym)
            SYNONYM_GROUP[synonym] = SYNONYM_GROUP[ab_name]
        else:
            SYNONYM_GROUP[ab_name] = set()
            SYNONYM_GROUP[ab_name].add(ab_name)
            SYNONYM_GROUP[ab_name].add(synonym)
            SYNONYM_GROUP[synonym] = SYNONYM_GROUP[ab_name]


AB_NAME2MAB_CLASS = {}


def init_abname2class(conn):
    cursor = conn.cursor()
    sql = """
        SELECT * FROM antibody_targets
    """

    cursor.execute(sql)

    for item in cursor.fetchall():
        ab_name = item[0]
        target = item[1]
        ab_class = item[2]
        source = item[3]
        class_info = {
            'target': target,
            'class': ab_class,
            'source': source
        }

        AB_NAME2MAB_CLASS[ab_name] = class_info

        for synonym in SYNONYM_GROUP.get(ab_name, []):
            AB_NAME2MAB_CLASS[synonym] = class_info


SUSCEPTIBLE_LEVEL_FILTER = """
    AND (
        (fold < 3
             OR (fold = 3 AND fold_cmp = '=')
             OR (fold = 3 AND fold_cmp = '<'))
        OR (resistance_level = 'susceptible'))
"""

PARTIAL_RESISTANCE_LEVEL_FILTER = """
    AND (
        (
            (fold = 3 AND fold_cmp = '>')
            OR (fold > 3 AND fold < 10)
            OR (fold = 10 AND fold_cmp = '=')
            OR (fold = 10 AND fold_cmp = '<'))
        OR (resistance_level = 'partial-resistance'))
"""

RESISTANCE_LEVLE_FILTER = """
    AND (
        (
            (fold = 10 AND fold_cmp = '>')
            OR (fold > 10))
        OR (resistance_level = 'partial-resistance'))
"""

RESISTANCE_FILTER = {
    'susceptible': [SUSCEPTIBLE_LEVEL_FILTER],
    'partial': [PARTIAL_RESISTANCE_LEVEL_FILTER],
    'resistance': [RESISTANCE_LEVLE_FILTER],
}
