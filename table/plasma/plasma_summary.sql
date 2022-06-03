SELECT
    *,
    julianday(collection_date) - julianday(last_event_date) days
FROM (
    SELECT
        plasma.ref_name,
        plasma.rx_name,
        plasma.subject_name,
        plasma.collection_date,
        hist.event_date,
        MAX(hist.event_date)
            OVER (PARTITION BY plasma.ref_name, plasma.rx_name) last_event_date,
        group_concat(hist.infection, "+")
            OVER (PARTITION BY plasma.ref_name, plasma.rx_name) infections,
        group_concat(hist.vaccine_name, "+")
            OVER (PARTITION BY plasma.ref_name, plasma.rx_name) vaccine_names,
        MAX(hist.dosage)
            OVER (PARTITION BY plasma.ref_name, plasma.rx_name) dosage
    FROM
        subject_plasma plasma,
        ({vp_plasma_history}) hist
    WHERE
        plasma.ref_name = hist.ref_name
        AND
        plasma.subject_name = hist.subject_name

        AND
        plasma.collection_date > hist.event_date
    )
WHERE
    last_event_date = event_date
;
