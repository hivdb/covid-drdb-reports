AGGREGATED_SUSC_VIEW_SQL = """
SELECT
    ref_name,
    rx_name,
    variant_name,
    control_variant_name,
    fold_cmp,
    fold
    cumulative_count
FROM susc_results
WHERE cumulative_count > 1
"""

INDIVIDUAL_SAMPLE_SQL = """
SELECT
 *
FROM
    susc_results AS a,
    (SELECT
        ref_name,
        rx_name,
        variant_name,
        control_variant_name
    FROM susc_results
    EXCEPT
    SELECT
        ref_name,
        rx_name,
        variant_name,
        control_variant_name
    FROM susc_results
    WHERE cumulative_count > 1
    ) AS b
ON
    a.ref_name = b.ref_name
    AND a.rx_name = b.rx_name
    AND a.variant_name = b.variant_name
    AND a.control_variant_name = b.control_variant_name
"""

AGGREGATED_SUSC_VIEW_SQL = """
SELECT
*
FROM susc_results
WHERE cumulative_count > 1
"""


AGGREGATED_SAMPLES = []


def get_aggregated_studies(conn):
    global AGGREGATED_SAMPLES

    cursor = conn.cursor()
    cursor.execute(AGGREGATED_SUSC_VIEW_SQL)

    result = []
    for row in cursor.fetchall():
        result.append(dict(row))

    AGGREGATED_SAMPLES = result

    return result


