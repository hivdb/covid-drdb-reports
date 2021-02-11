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
        OR (resistance_level = 'resistance'))
"""

RESISTANCE_FILTER = {
    'susceptible': [SUSCEPTIBLE_LEVEL_FILTER],
    'partial': [PARTIAL_RESISTANCE_LEVEL_FILTER],
    'resistance': [RESISTANCE_LEVLE_FILTER],
}


EXCLUDE_PLASMA = [
    'Hospitalized_Samples',
    'Mild_Samples',
]

PLASMA_RENAME = {
    'BNT162b2_3Weeks': 'BNT162b2_1M',
    'CP_5-33d': 'CP_1M',
    'CP_8mon': 'CP_8M',
    'Mod_1M': 'mRNA-1273_1M',
    'Moderna_36d': 'mRNA-1273_1M',
    'Moderna_D43': 'mRNA-1273_1M',
    'NVV_1M': 'NVX-CoV2373_1M',
    'Pfizer_BNT162b2_D28': 'BNT162b2_1M',
    'Pfizer-BioNTech': 'BNT162b2',
    'Moderna': 'mRNA-1273',
    'CP_13': 'CP',
    'CP_29': 'CP',
    'CP_35': 'CP',
    'CP_37': 'CP',
    'CP_BNT162b2': 'BNT162b2',
    'CP_ModerateIgG': 'CP',
    'CP_StrongIgG': 'CP',
    'CP_weakIgG': 'CP',
    'Subject_B_d26': 'CP',
    'Subject_C_d32': 'CP',
    'Subject_G_d18': 'CP',
    'Subject_I_d26': 'CP',
    'CP_BNT162b2_1M': 'BNT162b2_1M',
    'BNT162b2_2Weeks': 'BNT162b2_1M',
    'BNT162b2_4Weeks': 'BNT162b2_1M',
}

MAB_RENAME = {
    'LY-CoV555/CB6': 'Bamlanivimab/Etesevimab',
    'REGN10933/10987': 'Casirivimab/Imdevimab',
    'COV2-2196/2130': 'Cilgavimab/Tixagevimab',
}
