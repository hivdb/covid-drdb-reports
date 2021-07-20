AGGREGATED_SUSC_VIEW_SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.iso_name,
    s.control_iso_name,
    s.cumulative_count,
    s.fold_cmp,
    s.fold,
    s.resistance_level,
    s.ineffective
FROM
    susc_results AS s
WHERE
    s.potency_type IN ('IC50', 'NT50')
    AND
    (s.ref_name, s.rx_name, s.iso_name, s.control_iso_name) IN
    (
    SELECT
        ref_name,
        rx_name,
        iso_name,
        control_iso_name
    FROM
        susc_results
    WHERE
        cumulative_count > 1
    )
"""

INDIVIDUAL_SAMPLE_SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.iso_name,
    s.control_iso_name,
    s.cumulative_count,
    s.fold_cmp,
    s.fold,
    s.resistance_level,
    s.ineffective
FROM
    susc_results AS s
WHERE
    s.potency_type IN ('IC50', 'NT50')
    AND
    (s.ref_name, s.rx_name, s.iso_name, s.control_iso_name) IN
    (SELECT
        ref_name,
        rx_name,
        iso_name,
        control_iso_name
    FROM
        susc_results
    EXCEPT
    SELECT
        ref_name,
        rx_name,
        iso_name,
        control_iso_name
    FROM
        susc_results
    WHERE
        cumulative_count > 1
    )
"""
