SELECT
    iso.var_name,
    plasma.ref_name,
    plasma.rx_name,
    rx.potency titer,
    rx.cumulative_count num_exp,
    plasma.collection_date,
    plasma.event_date,
    plasma.days,
    plasma.infections,
    plasma.vaccine_names,
    plasma.dosage
FROM
    rx_potency rx,
    ({plasma_summary}) plasma,
    isolates iso,
    subjects subj
WHERE
    rx.ref_name = plasma.ref_name
    AND
    rx.rx_name = plasma.rx_name
    AND

    rx.potency_type = 'NT50'
    AND

    rx.iso_name = iso.iso_name
    AND
    iso.var_name IS NOT NULL
    AND

    plasma.ref_name = subj.ref_name
    AND
    plasma.subject_name = subj.subject_name
    AND
    subj.subject_species = 'Human'
;
