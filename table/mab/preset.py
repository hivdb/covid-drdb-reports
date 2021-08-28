ANTIBODY_TARGET_SQL = """
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
"""

ANTIBODY_EPITOPE_SQL = """
SELECT
    ab_name,
    GROUP_CONCAT(position, '+') AS epitope
FROM
    antibody_epitopes
GROUP BY
    ab_name
"""

ANTIBODY_SYNONYM_SQL = """
SELECT
    ab_name,
    GROUP_CONCAT(synonym, ';') AS synonyms
FROM
    antibody_synonyms
GROUP BY
    ab_name
"""

ANTIBODY_INFO_SQL = """
SELECT
    a.ab_name,
    s.synonyms,
    a.availability,
    t.pdb_id,
    t.target,
    t.class,
    e.epitope
FROM
    antibodies a
    LEFT JOIN
    ({antibody_target_sql}) t
    ON
        a.ab_name = t.ab_name
    LEFT JOIN
    ({antibody_epitope_sql}) e
    ON
        a.ab_name = e.ab_name
    LEFT JOIN
    ({antibody_synonym_sql}) s
    ON
        a.ab_name = s.ab_name
""".format(
    antibody_target_sql=ANTIBODY_TARGET_SQL,
    antibody_epitope_sql=ANTIBODY_EPITOPE_SQL,
    antibody_synonym_sql=ANTIBODY_SYNONYM_SQL)


RX_SINGLE_MAB_SQL = """
SELECT
    r.ref_name AS ref_name,
    r.rx_name AS rx_name,
    r.ab_name AS ab_name,
    ab.synonyms AS synonyms,
    ab.availability,
    ab.pdb_id,
    ab.target,
    ab.class,
    ab.epitope
FROM
    rx_antibodies AS r,
    ({antibody_info_sql}) AS ab
ON r.ab_name = ab.ab_name
GROUP BY r.ref_name, r.rx_name
HAVING count(r.ab_name) = 1

UNION

SELECT
    r.ref_name AS ref_name,
    r.rx_name AS rx_name,
    r.ab_name AS ab_name,
    ab.synonyms AS synonyms,
    ab.availability,
    ab.pdb_id,
    ab.target,
    ab.class,
    ab.epitope
FROM
    rx_dms AS r,
    ({antibody_info_sql}) AS ab
ON r.ab_name = ab.ab_name
GROUP BY r.ref_name, r.rx_name
HAVING count(r.ab_name) = 1

""".format(
    antibody_info_sql=ANTIBODY_INFO_SQL)

RX_COMBO_MAB_SQL = """
SELECT
    ref_name,
    rx_name,
    group_concat(a.ab_name, '/') AS ab_name,
    '' AS synonyms,
    availability AS availability,
    '' AS pdb_id,
    '' AS target,
    '' AS class,
    '' AS epitope
FROM
    (SELECT * FROM rx_antibodies ORDER BY ref_name, ab_name) AS a,
    antibodies AS b
ON
    a.ab_name = b.ab_name
GROUP BY
    ref_name, rx_name
HAVING
    count(a.ab_name) > 1

UNION

SELECT
    ref_name,
    rx_name,
    group_concat(a.ab_name, '/') AS ab_name,
    '' AS synonyms,
    availability AS availability,
    '' AS pdb_id,
    '' AS target,
    '' AS class,
    '' AS epitope
FROM
    (SELECT * FROM rx_dms ORDER BY ref_name, ab_name) AS a,
    antibodies AS b
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
print(RX_MAB)


MAB_RENAME = {}


def load_mab_rename(conn):

    cursor = conn.cursor()

    cursor.execute(RX_MAB)

    for row in cursor.fetchall():
        rx_name = row['rx_name']
        ab_name = row['ab_name']
        MAB_RENAME[rx_name] = ab_name
