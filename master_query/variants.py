def in_clause(in_list):
    return ', '.join(["'{}'".format(i) for i in in_list])


WT_VARIANTS = """
SELECT
    iso_name
FROM
    isolates
WHERE
    var_name in ('A', 'B', 'B.1')
    OR iso_name in ('S:614G', 'S:683G')
"""

IGNORE_ISOLATES = [
    'SARS-CoV Spike',
    'WIV1 Spike',
]

SPIKE_ONLY_ISOLATES = """
SELECT
    *
FROM
    isolate_mutations
WHERE
    iso_name NOT IN
    (
        SELECT
            DISTINCT iso_name
        FROM
            isolate_mutations
        WHERE
            gene != 'S'
    )
    AND
    iso_name NOT IN ({exclude_iso})
""".format(
    exclude_iso=in_clause(IGNORE_ISOLATES)
)

SINGLE_SPIKE_MUTATION_ISOLATES = """
SELECT
    iso_name
FROM
    ({spike_only_isolates})
GROUP BY
    iso_name
HAVING
    count(position) = 1

UNION

SELECT
    iso_name
FROM
    ({spike_only_isolates})
GROUP BY
    iso_name
HAVING
    count(position) = 2
    AND
    (
        group_concat(position) LIKE '%614%'
        OR
        group_concat(position) LIKE '%683%'
    )
""".format(
    spike_only_isolates=SPIKE_ONLY_ISOLATES
)


SINGLE_SPIKE_MUTATION = """
SELECT
    iso_name,
    position,
    amino_acid
FROM
    isolate_mutations
WHERE
    iso_name in ({single_spike_mutation_isolates})
    AND
    position not in (614, 683)
""".format(
    single_spike_mutation_isolates=SINGLE_SPIKE_MUTATION_ISOLATES
    )
