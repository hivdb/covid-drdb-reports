INDIVIDUAL_CP_SQL = """
SELECT
    s.ref_name as ref_name,
    rx.cumulative_group as rx_name,
    rx.timing,
    rx.severity,
    s.control_iso_name as control,
    s.iso_name as iso_name,
    s.cumulative_count as num_fold,
    s.fold_cmp as fold_cmp,
    s.fold as fold
FROM
    susc_results_indiv_wt_view s,
    {rx_type} AS rx
WHERE
    rx.ref_name = s.ref_name
    AND rx.rx_name = s.rx_name
"""


INDIVIDUAL_VP_SQL = """
SELECT
    s.ref_name as ref_name,
    rx.vaccine_name as rx_name,
    rx.timing,
    rx.dosage,
    s.control_iso_name as control,
    s.iso_name as iso_name,
    s.cumulative_count as num_fold,
    s.fold_cmp as fold_cmp,
    s.fold as fold
FROM
    susc_results_indiv_wt_view s,
    {rx_type} AS rx
WHERE
    rx.ref_name = s.ref_name
    AND rx.rx_name = s.rx_name
"""


AGGREGATED_CP_SQL = """
SELECT
    s.ref_name as ref_name,
    rx.cumulative_group as rx_name,
    rx.timing,
    rx.severity,
    s.control_iso_name as control,
    s.iso_name as iso_name,
    s.cumulative_count as num_fold,
    s.fold_cmp as fold_cmp,
    s.fold as fold
FROM
    susc_results_aggr_wt_view s,
    {rx_type} AS rx
WHERE
    rx.ref_name = s.ref_name
    AND rx.rx_name = s.rx_name
"""

AGGREGATED_VP_SQL = """
SELECT
    s.ref_name as ref_name,
    rx.vaccine_name as rx_name,
    rx.timing,
    rx.dosage,
    s.control_iso_name as control,
    s.iso_name as iso_name,
    s.cumulative_count as num_fold,
    s.fold_cmp as fold_cmp,
    s.fold as fold
FROM
    susc_results_aggr_wt_view s,
    {rx_type} AS rx
WHERE
    rx.ref_name = s.ref_name
    AND rx.rx_name = s.rx_name
"""


IGNORE_VACCINE_NAME = [
    'mRNA-1273.211',
    'mRNA-1273.351',
    'mRNA-1273 + mRNA-1273.351',
]

IGNORE_STUDY = [
    'Wall21b'
]


INFECTED_VACCINEE = """
SELECT DISTINCT ref_name,
                subject_name
FROM   subject_history a
WHERE  EXISTS (SELECT 1
               FROM   subject_history b
               WHERE  a.subject_name = b.subject_name
                      AND a.ref_name = b.ref_name
                      AND b.event = 'infection')
       AND EXISTS (SELECT 1
                   FROM   subject_history b
                   WHERE  a.subject_name = b.subject_name
                          AND a.ref_name = b.ref_name
                          AND b.event LIKE '%dose%')
"""


RX_TITER_FOLD = """
SELECT *
FROM susc_results a,
     rx_potency b
WHERE
    a.ref_name = b.ref_name
  AND a.rx_name = b.rx_name
  AND a.iso_name = b.iso_name
"""
