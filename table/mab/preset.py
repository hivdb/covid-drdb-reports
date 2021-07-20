MAB_RENAME = {
    'LY-CoV555/CB6': 'Bamlanivimab/Etesevimab',
    'LY-CoV555+LY-CoV016': 'Bamlanivimab/Etesevimab',
    'Bamlanivimab+Etesevimab': 'Bamlanivimab/Etesevimab',

    'COV2-2196/2130': 'Cilgavimab/Tixagevimab',
    'Tixagevimab + Cilgavimab': 'Cilgavimab/Tixagevimab',
    'AZD7442': 'Cilgavimab/Tixagevimab',
    'COV2-2130+COV2-2196': 'Cilgavimab/Tixagevimab',

    'REGN10933/10987': 'Casirivimab/Imdevimab',
    'REGN10933+REGN10987': 'Casirivimab/Imdevimab',
    'REGN10933 + REGN10987': 'Casirivimab/Imdevimab',
    'REGN10933+REGN10987': 'Casirivimab/Imdevimab',
    'Casirivimab+Imdevimab': 'Casirivimab/Imdevimab',
    'Casirivimab + Imdevimab': 'Casirivimab/Imdevimab',
    'CAS/IMD': 'Casirivimab/Imdevimab',

    'BRII-196/198': 'BRII-196/BRII-198',

    'C144+C135': 'C135/C144',

    'Vir-7831+S2E12': 'Vir-7831/S2E12',

    'REGN10989/10987': 'REGN10989/Imdevimab',
    'REGN10989+10934': 'REGN10989+10934',

    'LY-CoV555': 'Bamlanivimab',
    'LY3819253': 'Bamlanivimab',
    'REGN10933': 'Casirivimab',
    'CB6': 'Etesevimab',
    'JS-016': 'Etesevimab',
    'LY-CoV016': 'Etesevimab',
    'LY3832479': 'Etesevimab',
    'REGN10987': 'Imdevimab',
    'GSK4182136': 'Vir-7831',
    'S309': 'Vir-7831',
    'Sotrovimab': 'Vir-7831',
    'COV2-2130': 'Cilgavimab',
    'AZD1061': 'Cilgavimab',
    'COV2-2196': 'Tixagevimab',
    'AZD8895': 'Tixagevimab',
    'CT-P59': 'Regdanvimab',

}

ANTIBODY_TARGET_SQL = """
(
    SELECT
        DISTINCT *
    FROM
        antibody_targets
    WHERE
        ab_name NOT IN (
            SELECT
                ab_name
            FROM
                antibody_targets
            WHERE
                pdb_id IS NOT null
            )
UNION
    SELECT
        DISTINCT *
    FROM
        antibody_targets
    WHERE
        pdb_id IS NOT null
)
"""

ANTIBODY_INFO_SQL = """
SELECT
    a.ab_name,
    a.availability,
    b.pdb_id,
    b.target,
    b.class
FROM
    antibodies AS a
    LEFT JOIN
    {antibody_target_sql} AS b
ON
    a.ab_name = b.ab_name
""".format(antibody_target_sql=ANTIBODY_TARGET_SQL)

RX_SINGLE_MAB_SQL = """
SELECT
    r.ref_name as ref_name,
    r.rx_name as rx_name,
    r.ab_name as ab_name,
    ab.availability,
    ab.pdb_id,
    ab.target,
    ab.class
FROM
    rx_antibodies as r,
    ({antibody_info_sql}) as ab
ON r.ab_name = ab.ab_name
GROUP BY r.ref_name, r.rx_name
HAVING count(r.ab_name) = 1
""".format(antibody_info_sql=ANTIBODY_INFO_SQL)

RX_COMBO_MAB_SQL = """
SELECT
    ref_name,
    rx_name,
    group_concat(a.ab_name, '/') as ab_name,
    availability as availability,
    '' as pdb_id,
    '' as target,
    '' as class
FROM
    (SELECT * FROM rx_antibodies ORDER BY ref_name, ab_name) as a,
    antibodies as b
ON
    a.ab_name = b.ab_name
GROUP BY
    ref_name, rx_name
HAVING
    count(a.ab_name) > 1
"""

RX_MAB = """{single_mab} UNION {combo_mab}""".format(
    single_mab=RX_SINGLE_MAB_SQL,
    combo_mab=RX_COMBO_MAB_SQL)
