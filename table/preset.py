import csv
from pathlib import Path
import decimal
from decimal import Decimal
import json


decimal.getcontext().rounding = decimal.ROUND_HALF_UP

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


def load_csv(file_path):
    records = []
    with open(file_path, encoding='utf-8-sig') as fd:
        for record in csv.DictReader(fd):
            records.append(record)

    return records


def dump_json(json_path, obj):
    with json_path.open('w') as fd:
        json.dump(obj, fd, indent=4, ensure_ascii=False)


SYNONYM2AB_NAME = {}
AB_NAME2SYNONYM = {}
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
        AB_NAME2SYNONYM[ab_name] = synonym

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

    cursor = conn.cursor()
    sql = """
        SELECT ab_name, abbreviation_name
        FROM antibodies WHERE abbreviation_name IS NOT NULL;
    """
    cursor.execute(sql)

    for item in cursor.fetchall():
        ab_name = item[0]
        abbr_name = item[1]
        if ab_name in AB_NAME2SYNONYM.keys():
            SYNONYM2AB_NAME[abbr_name] = ab_name
        elif ab_name in SYNONYM2AB_NAME.keys():
            SYNONYM2AB_NAME[abbr_name] = SYNONYM2AB_NAME[ab_name]

        SYNONYM_GROUP[ab_name].add(abbr_name)


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
             OR (fold = 3 AND fold_cmp = '<')
             OR (fold = 3 AND fold_cmp = '=')
             OR (fold = 3 AND fold_cmp = '~'))
        OR (resistance_level = 'susceptible')
        )
    AND (ineffective IS NULL)
"""

PARTIAL_RESISTANCE_LEVEL_FILTER = """
    AND (
        (
            (fold = 3 AND fold_cmp = '>')
            OR (fold > 3 AND fold < 10)
            OR (fold = 10 AND fold_cmp = '=')
            OR (fold = 10 AND fold_cmp = '~')
            OR (fold = 10 AND fold_cmp = '<'))
        OR (resistance_level = 'partial-resistance'))
    AND (ineffective IS NULL)
"""

RESISTANT_LEVLE_FILTER = """
    AND (
        (
            (fold = 10 AND fold_cmp = '>')
            OR (fold > 10))
        OR (resistance_level = 'resistant')
        OR (ineffective = 'experimental')
        )
"""

RESISTANCE_FILTER = {
    'susceptible': [SUSCEPTIBLE_LEVEL_FILTER],
    'partial': [PARTIAL_RESISTANCE_LEVEL_FILTER],
    'resistant': [RESISTANT_LEVLE_FILTER],
}


EXCLUDE_STUDIES = {
    'Planas21': lambda x: x.startswith('CP') or (not x.endswith('W4')),
    'Garcia-Beltran21': lambda x: x.startswith('BNT') or x.startswith('mRNA'),
    'Collier21': lambda x: x.startswith('CP'),
    'Widera21': lambda x: x.startswith('CP') or x.startswith('BNT'),
    'Brown21': lambda x: x.startswith('CP') or x.startswith('BNT'),
    'Supasa21': lambda x: x.startswith('CP')
                          or x.startswith('BNT') or x.startswith('AZD'),
    'Edara21': lambda x: x == 'CP_acute',
}


EXCLUDE_PLASMA = [
    'Severe',
    'Mild',
]


RENAME_CP_EXECUTOR = {
    'Planas21': [
        (
            lambda x: x.startswith('VAC') and x.endswith('W4'),
            'BNT162b2'
        ),
        (
            lambda x: x.startswith('VAC') and x.endswith('W2'),
            'BNT162b2_W2'
        ),
        (
            lambda x: x.startswith('VAC') and x.endswith('W3'),
            'BNT162b2_W3'
        ),
    ],
    'Hoffmann21b': [
        (
            lambda x: x.startswith('ConvSerum') or x.startswith('ConvPlasma'),
            'CP',
        ),
    ]
}


PLASMA_RENAME = {
    'BNT162b2_3W': 'BNT162b2_1M',
    'CP_5-33d': 'CP_1M',
    'CP_8M': 'CP_8M',
    'Mod_1M': 'mRNA-1273_1M',
    'Moderna_36d': 'mRNA-1273_1M',
    'Moderna_D43': 'mRNA-1273_1M',
    'NVV_1M': 'NVX-CoV2373_1M',
    'Pfizer_BNT162b2_D28': 'BNT162b2_1M',
    'Pfizer-BioNTech': 'BNT162b2',
    'Moderna': 'mRNA-1273',
    'BNT162b2_infected': 'BNT162b2',
    'CP_Mild': 'CP',
    "CP_Patient": 'CP',
    'CP_13': 'CP',
    'CP_29': 'CP',
    'CP_35': 'CP',
    'CP_37': 'CP',
    'CP_614G': 'CP',
    'CP_BNT162b2': 'BNT162b2',
    'CP_ModerateIgG': 'CP',
    'CP_StrongIgG': 'CP',
    'CP_weakIgG': 'CP',
    'Subject_B_d26': 'CP',
    'Subject_C_d32': 'CP',
    'Subject_G_d18': 'CP',
    'Subject_I_d26': 'CP',
    'CP_BNT162b2_1M': 'BNT162b2_1M',
    'BNT162b2_2W': 'BNT162b2_1M',
    'BNT162b2_4W': 'BNT162b2_1M',
    'CP_02-0014': 'CP_WT',
    'CP_02-0015': 'CP_WT',
    'CP_13-0013': 'CP_WT',
    'CP_13-0017': 'CP_WT',
    'CP_13-0033': 'CP_WT',
    'CP_13-0062': 'CP_WT',
    'BNT162b2_1M:01': 'BNT162b2',
    'BNT162b2_1M:02': 'BNT162b2',
    'BNT162b2_1M:03': 'BNT162b2',
    'BNT162b2_1M:04': 'BNT162b2',
    'BNT162b2_1M:05': 'BNT162b2',
    'BNT162b2_1M:06': 'BNT162b2',
    'BNT162b2_1M:07': 'BNT162b2',
    'BNT162b2_1M:08': 'BNT162b2',
    'BNT162b2_1M:09': 'BNT162b2',
    'BNT162b2_1M:10': 'BNT162b2',
    'BNT162b2_1M:11': 'BNT162b2',
    'BNT162b2_1M:12': 'BNT162b2',
    'BNT162b2_1M:13': 'BNT162b2',
    'BNT162b2_1M:14': 'BNT162b2',
    'BNT162b2_1M:15': 'BNT162b2',
    'CP_ID15': 'CP',
    'CP_ID18': 'CP',
    'CP_ID20': 'CP',
    'CP_ID22': 'CP',
    'CP_ID23': 'CP',
    'CP_ID24': 'CP',
    'CP_ID27': 'CP',
    'CP_ID33': 'CP',
    'CP_ID51': 'CP',
    'hCoV-2IG': 'IVIG',
}

PLASMA_POST_RENAME = {
    'BNT162b2_1M': 'BNT162b2',
    'CP_1M': 'CP',
    'CP_8M': 'CP(8M)',
    'CP_6M': 'CP(6M)',
    'NVX-CoV2373_1M': 'NVX-CoV',
    'mRNA-1273_1M': 'mRNA-1273',
    'CP_Titer_GT400': 'CP(High Titer)',
    'CP_Titer_LT400': 'CP(Low Titer)',
    'CP_WT': 'CP(WT)',
}

MAB_RENAME = {
    'LY-CoV555/CB6': 'Bamlanivimab/Etesevimab',
    'REGN10933/10987': 'Casirivimab/Imdevimab',
    'COV2-2196/2130': 'Cilgavimab/Tixagevimab',
    'REGN10933 + REGN10987': 'Casirivimab/Imdevimab',
    'Vir-7831': 'Sotrovimab',
}

EXCLUDE_MAB = [
    'BRII-196',
    'BRII-198',
    'REGN10989/10987',
    'BRII-196/198',
    'DH1041',
    '5-24',
    '910-30',
    'BRII-196+BRII-198',
]


def round_number(number):
    if number < 10:
        number = float(number)
        if number.is_integer():
            return Decimal(str(number)).quantize(Decimal('1'))
        else:
            return Decimal(str(number)).quantize(Decimal('1.0'))
    elif number >= 10 and number < 100:
        return Decimal(str(number)).quantize(Decimal('1'))
    else:
        return '>100'
