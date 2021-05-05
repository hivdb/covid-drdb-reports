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
    susc_results AS s
WHERE
    s.inhibition_pcnt != 90
    AND
    (s.ref_name, s.rx_name, s.variant_name, s.control_variant_name) IN
    (
    SELECT
        ref_name,
        rx_name,
        variant_name,
        control_variant_name
    FROM susc_results
    WHERE cumulative_count > 1
    )
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
    susc_results AS s
WHERE
    s.inhibition_pcnt != 90
    AND
    (s.ref_name, s.rx_name, s.variant_name, s.control_variant_name) IN
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
    )
"""
