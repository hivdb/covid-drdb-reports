SELECT * FROM (
    SELECT
        ref_name,
        subject_name,
        date(julianday(infection_date) - 1) event_date,  -- test positive and get sample the same day.
        infected_var_name infection,
        NULL vaccine_name,
        NULL dosage
    FROM
        subject_infections

    UNION ALL

    SELECT
        ref_name,
        subject_name,
        vaccination_date event_date,
        NULL infection,
        vaccine_name,
        dosage
    FROM
        subject_vaccines
)
ORDER BY
ref_name,
subject_name,
event_date asc
;
