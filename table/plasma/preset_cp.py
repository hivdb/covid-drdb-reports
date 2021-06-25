CP_VIEW_INFECTION = """
SELECT
    a.ref_name,
    a.rx_name,
    event_date as infection_date,
    iso.var_name,
    c.severity
FROM
    rx_conv_plasma AS a,
    patient_treatments AS b,
    patient_history AS c,
    isolates AS iso
ON
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    b.ref_name = c.ref_name
    AND
    b.patient_name = c.patient_name
    AND
    c.iso_name = iso.iso_name
WHERE
    c.event = 'infection'
"""

CP_VIEW_ISOLATION = """
SELECT
    a.ref_name,
    a.rx_name,
    a.collection_date as isolation_date,
    iso.var_name,
    c.severity
FROM
    rx_conv_plasma AS a,
    patient_treatments AS b,
    patient_history AS c,
    isolates AS iso
ON
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    b.ref_name = c.ref_name
    AND
    b.patient_name = c.patient_name
    AND
    c.iso_name = iso.iso_name
WHERE
    c.event = 'isolation'
    AND
    c.event_date = a.collection_date
"""

CP_VIEW_WITH_HISTORY = """
SELECT
    a.ref_name,
    a.rx_name,
    b.infection_date,
    a.isolation_date,
    ROUND(
        (
            julianday(a.isolation_date) -
            julianday(b.infection_date)
        ) / 30, 1) as timing,
    a.var_name,
    a.severity
FROM
    ({isolation_view}) AS a,
    ({infection_view}) AS b
ON
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    a.var_name = b.var_name
    AND
    a.severity = b.severity
UNION
SELECT
    a.ref_name,
    a.rx_name,
    b.infection_date,
    a.isolation_date,
    NULL as timing,
    a.var_name,
    a.severity
FROM
    ({isolation_view}) AS a
    LEFT JOIN
    ({infection_view}) AS b
ON
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    a.var_name = b.var_name
    AND
    a.severity = b.severity
WHERE
    b.infection_date IS NULL
""".format(
    infection_view=CP_VIEW_INFECTION,
    isolation_view=CP_VIEW_ISOLATION
)

CP_VIEW_WITHOUT_HISTORY = """
SELECT
    a.ref_name,
    a.rx_name,
    NULL AS infection_date,
    NULL AS isolation_date,
    NULL AS timing,
    'Unknown' AS var_name,
    NULL AS severity
FROM
    rx_conv_plasma AS a
    LEFT JOIN
    patient_treatments AS b
ON
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
WHERE
    b.patient_name IS NULL
"""

CP_VIEW = """
{with_history}
UNION
{without_hisotry}
""".format(
    with_history=CP_VIEW_WITH_HISTORY,
    without_hisotry=CP_VIEW_WITHOUT_HISTORY
)
