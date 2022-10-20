SELECT
    s.ref_name,
    rx.ab_name,
    rx.class,
    s.fold_cmp,
    s.fold,
    s.ineffective
FROM
    {susc_table} as s,
    rx_mab_view as rx,
    isolate_mutations_combo_s_mut_view mut
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.iso_name = mut.iso_name
    AND
    s.fold IS NOT NULL
    AND (
        rx.availability IS NOT NULL
        OR rx.pdb_id IS NOT NULL
    )
    AND
    {filters}
