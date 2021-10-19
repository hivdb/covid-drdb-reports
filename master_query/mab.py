MAB_TARGET_SQL = """
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

MAB_EPITOPE_SQL = """
SELECT
    ab_name,
    GROUP_CONCAT(position, '+') AS epitope
FROM
    antibody_epitopes
GROUP BY
    ab_name
"""

MAB_SYNONYM_SQL = """
SELECT
    ab_name,
    GROUP_CONCAT(synonym, ';') AS synonyms
FROM
    antibody_synonyms
GROUP BY
    ab_name
"""

SINGLE_MAB_INFO_SQL = """
SELECT
    a.ab_name,
    s.synonyms,
    a.availability,
    t.pdb_id,
    t.target,
    t.class,
    e.epitope,
    a.institute,
    a.origin
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
    antibody_target_sql=MAB_TARGET_SQL,
    antibody_epitope_sql=MAB_EPITOPE_SQL,
    antibody_synonym_sql=MAB_SYNONYM_SQL)

RX_SINGLE_MAB_DRDB_SQL = """
SELECT
    r.ref_name AS ref_name,
    r.rx_name AS rx_name,
    r.ab_name AS ab_name,
    ab.synonyms AS synonyms,
    ab.availability,
    ab.pdb_id,
    ab.target,
    ab.class,
    ab.epitope,
    ab.institute,
    ab.origin
FROM
    rx_antibodies AS r,
    ({single_mab_info_sql}) AS ab
ON r.ab_name = ab.ab_name
GROUP BY r.ref_name, r.rx_name
HAVING count(r.ab_name) = 1
""".format(
    single_mab_info_sql=SINGLE_MAB_INFO_SQL
    )

RX_SINGLE_MAB_DMS_SQL = """
SELECT
    r.ref_name AS ref_name,
    r.rx_name AS rx_name,
    r.ab_name AS ab_name,
    ab.synonyms AS synonyms,
    ab.availability,
    ab.pdb_id,
    ab.target,
    ab.class,
    ab.epitope,
    ab.institute,
    ab.origin
FROM
    rx_antibodies r,
    ({single_mab_info_sql}) AS ab
ON r.ab_name = ab.ab_name
GROUP BY r.ref_name, r.rx_name
HAVING count(r.ab_name) = 1
""".format(
    single_mab_info_sql=SINGLE_MAB_INFO_SQL
    )

RX_SINGLE_MAB_SQL = """
{rx_single_mab_fold}
UNION
{rx_single_mab_dms}
""".format(
    rx_single_mab_fold=RX_SINGLE_MAB_DRDB_SQL,
    rx_single_mab_dms=RX_SINGLE_MAB_DMS_SQL
    )

RX_COMBO_MAB_DRDB_SQL = """
SELECT
    a.ref_name,
    a.rx_name,
    a.ab_name || '/' || b.ab_name AS ab_name,
    NULL AS synonyms,
    CASE
        WHEN a.availability == b.availability THEN
            a.availability
        ELSE
            NULL
    END availability,
    NULL AS pdb_id,
    NULL AS target,
    NULL AS class,
    NULL AS epitope,
    CASE
        WHEN a.institute == b.institute THEN
            a.institute
        ELSE
            NULL
    END institute,
    CASE
        WHEN a.origin == b.origin THEN
            a.origin
        ELSE
            NULL
    END origin
FROM
    (SELECT
        a.*,
        b.availability,
        b.institute,
        b.origin
    FROM
        rx_antibodies a,
        antibodies b
    ON
        a.ab_name = b.ab_name
    ) a,
    (SELECT
        a.*,
        b.availability,
        b.institute,
        b.origin
    FROM
        rx_antibodies a,
        antibodies b
    ON
        a.ab_name = b.ab_name
    ) b
WHERE
    a.ref_name = b.ref_name and
    a.rx_name = b.rx_name and
    a.ab_name != b.ab_name and
    a.ab_name < b.ab_name
"""


RX_COMBO_MAB_DMS_SQL = """
SELECT
    a.ref_name,
    a.rx_name,
    a.ab_name || '/' || b.ab_name AS ab_name,
    NULL AS synonyms,
    CASE
        WHEN a.availability == b.availability THEN
            a.availability
        ELSE
            NULL
    END availability,
    NULL AS pdb_id,
    NULL AS target,
    NULL AS class,
    NULL AS epitope,
    CASE
        WHEN a.institute == b.institute THEN
            a.institute
        ELSE
            NULL
    END institute,
    CASE
        WHEN a.origin == b.origin THEN
            a.origin
        ELSE
            NULL
    END origin
FROM
    (SELECT
        a.*,
        b.availability,
        b.institute,
        b.origin
    FROM
        rx_antibodies a,
        antibodies b
    ON
        a.ab_name = b.ab_name
    ) a,
    (SELECT
        a.*,
        b.availability,
        b.institute,
        b.origin
    FROM
        rx_antibodies a,
        antibodies b
    ON
        a.ab_name = b.ab_name
    ) b
WHERE
    a.ref_name = b.ref_name and
    a.rx_name = b.rx_name and
    a.ab_name != b.ab_name and
    a.ab_name < b.ab_name
"""

RX_COMBO_MAB_SQL = """
{rx_combo_mab_fold}

UNION

{rx_combo_mab_dms}
""".format(
    rx_combo_mab_fold=RX_COMBO_MAB_DRDB_SQL,
    rx_combo_mab_dms=RX_COMBO_MAB_DMS_SQL
)

RX_MAB = """
{single_mab}
UNION
{combo_mab}
""".format(
    single_mab=RX_SINGLE_MAB_SQL,
    combo_mab=RX_COMBO_MAB_SQL
    )

RX_MAB_DRDB = """
{single_mab}
UNION
{combo_mab}
""".format(
    single_mab=RX_SINGLE_MAB_DRDB_SQL,
    combo_mab=RX_COMBO_MAB_DRDB_SQL
    )


RX_MAB_DMS = """
{rx_single_mab_dms}
UNION
{rx_combo_mab_dms}
""".format(
    rx_single_mab_dms=RX_SINGLE_MAB_DMS_SQL,
    rx_combo_mab_dms=RX_COMBO_MAB_DMS_SQL
    )
