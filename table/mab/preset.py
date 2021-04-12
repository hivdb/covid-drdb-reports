from pprint import pprint

SYNONYM2AB_NAME = {}
AB_NAME2SYNONYM = {}
SYNONYM_GROUP = {}


def init_synonyms_map(conn):
    cursor = conn.cursor()
    sql = """
        SELECT * FROM antibody_synonyms
    """
    cursor.execute(sql)

    for row in cursor.fetchall():
        ab_name = row['ab_name']
        synonym = row['synonym']
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

    for row in cursor.fetchall():
        ab_name = row['ab_name']
        abbr_name = row['abbreviation_name']
        if ab_name in AB_NAME2SYNONYM.keys():
            SYNONYM2AB_NAME[abbr_name] = ab_name
        elif ab_name in SYNONYM2AB_NAME.keys():
            SYNONYM2AB_NAME[abbr_name] = SYNONYM2AB_NAME[ab_name]

        SYNONYM_GROUP[ab_name].add(abbr_name)

    MAB_RENAME.update(SYNONYM2AB_NAME)


AB_NAME2MAB_CLASS = {}


def init_abname2class(conn):
    cursor = conn.cursor()
    sql = """
        SELECT * FROM antibody_targets
    """

    cursor.execute(sql)

    for row in cursor.fetchall():
        ab_name = row['ab_name']
        target = row['target']
        ab_class = row['class']
        source = row['source']
        class_info = {
            'target': target,
            'class': ab_class,
            'source': source
        }

        AB_NAME2MAB_CLASS[ab_name] = class_info

        for synonym in SYNONYM_GROUP.get(ab_name, []):
            AB_NAME2MAB_CLASS[synonym] = class_info


MAB_RENAME = {
    'LY-CoV555/CB6': 'Bamlanivimab/Etesevimab',
    'LY-CoV555+LY-CoV016': 'Bamlanivimab/Etesevimab',
    'Bamlanivimab+Etesevimab': 'Bamlanivimab/Etesevimab',
    'COV2-2196/2130': 'Cilgavimab/Tixagevimab',
    'COV2-2196+COV2-2130': 'Cilgavimab/Tixagevimab',
    'COV2-2130+COV2-2196': 'Cilgavimab/Tixagevimab',
    'Tixagevimab + Cilgavimab': 'Cilgavimab/Tixagevimab',
    'REGN10933 + REGN10987': 'Casirivimab/Imdevimab',
    'Casirivimab+Imdevimab': 'Casirivimab/Imdevimab',
    'Casirivimab + Imdevimab': 'Casirivimab/Imdevimab',
    'REGN10933/10987': 'Casirivimab/Imdevimab',
    'REGN10933+REGN10987': 'Casirivimab/Imdevimab',
    'BRII-196+BRII-198': 'BRII-196/BRII-198',
    'BRII-196/198': 'BRII-196/BRII-198',
}

EXCLUDE_MAB = [
    'REGN10989/10987',
    'DH1041',
    '5-24',
    '910-30',
]


ANTIBODY_TARGET_SQL = """
(SELECT DISTINCT * FROM antibody_targets
WHERE ab_name NOT IN (
    SELECT ab_name FROM 'antibody_targets' WHERE pdb_id IS NOT null)
UNION
SELECT DISTINCT * FROM antibody_targets
WHERE pdb_id IS NOT null)
"""
