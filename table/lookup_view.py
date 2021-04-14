AGGREGATED_SUSC_VIEW_SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.variant_name,
    s.control_variant_name,
    s.cumulative_count,
    s.fold_cmp,
    s.fold,
    s.resistance_level,
    s.ineffective
FROM
    susc_results AS s,
    (
    SELECT
        ref_name,
        rx_name,
        variant_name,
        control_variant_name
    FROM susc_results
    WHERE cumulative_count > 1
    ) AS b
ON
    s.ref_name = b.ref_name
    AND s.rx_name = b.rx_name
    AND s.variant_name = b.variant_name
    AND s.control_variant_name = b.control_variant_name
"""

INDIVIDUAL_SAMPLE_SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.variant_name,
    s.control_variant_name,
    s.cumulative_count,
    s.fold_cmp,
    s.fold,
    s.resistance_level,
    s.ineffective
FROM
    susc_results AS s,
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
    s.ref_name = b.ref_name
    AND s.rx_name = b.rx_name
    AND s.variant_name = b.variant_name
    AND s.control_variant_name = b.control_variant_name
"""
