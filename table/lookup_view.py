AGGREGATED_STUDY_SQL = """
SELECT
    ref_name,
    rx_name,
    variant_name,
    control_variant_name,
    cumulative_count
FROM susc_results
WHERE cumulative_count > 1
"""


def get_aggregated_studies(conn):
    cursor = conn.cursor()
    cursor.execute(AGGREGATED_STUDY_SQL)

    result = []
    for row in cursor.fetchall():
        result.append(dict(row))

    return result


