-- #Mab susc results
SELECT COUNT(*) FROM
    susc_results AS s, rx_antibodies AS rxmab
    WHERE rxmab.ref_name = s.ref_name AND rxmab.rx_name = s.rx_name;
-- #Mab susc results for one mutation
SELECT COUNT(*) FROM
    susc_results AS s,
    rx_antibodies AS rxmab,
    virus_strains AS vs,
    (SELECT strain_name, COUNT(*) AS num_muts FROM strain_mutations GROUP BY strain_name) AS sm
    WHERE rxmab.ref_name = s.ref_name AND rxmab.rx_name = s.rx_name
    AND vs.strain_name = s.strain_name
    AND sm.num_muts = 1 AND sm.strain_name = s.strain_name
    AND vs.site_directed IS TRUE;
-- #Mab susc results for mutation combinations (not strain)
SELECT COUNT(*) FROM
    susc_results AS s,
    rx_antibodies AS rxmab,
    virus_strains AS vs,
    (SELECT strain_name, COUNT(*) AS num_muts FROM strain_mutations GROUP BY strain_name) AS sm
    WHERE rxmab.ref_name = s.ref_name AND rxmab.rx_name = s.rx_name
    AND vs.strain_name = s.strain_name
    AND sm.num_muts > 1 AND sm.strain_name = s.strain_name
    AND vs.site_directed IS TRUE;
-- #Mab susc result for strains(not SDM)
SELECT COUNT(*) FROM
    susc_results AS s,
    rx_antibodies AS rxmab,
    virus_strains AS vs,
    (SELECT strain_name, COUNT(*) AS num_muts FROM strain_mutations GROUP BY strain_name) AS sm
    WHERE rxmab.ref_name = s.ref_name AND rxmab.rx_name = s.rx_name
    AND vs.strain_name = s.strain_name
    AND sm.strain_name = s.strain_name
    AND vs.site_directed IS FALSE
-- #CP susc results
SELECT COUNT(*) FROM susc_results AS s, rx_conv_plasma AS rxcp
    WHERE rxcp.ref_name = s.ref_name AND rxcp.rx_name = s.rx_name;
-- #CP susc results for one mutation
SELECT COUNT(*) FROM
    susc_results AS s,
    rx_conv_plasma AS rxcp,
    virus_strains AS vs,
    (SELECT strain_name, COUNT(*) AS num_muts FROM strain_mutations GROUP BY strain_name) AS sm
    WHERE rxcp.ref_name = s.ref_name AND rxcp.rx_name = s.rx_name
    AND vs.strain_name = s.strain_name
    AND sm.num_muts = 1 AND sm.strain_name = s.strain_name
    AND vs.site_directed IS TRUE
-- #CP susc results for mutation combinations(not strain)
SELECT COUNT(*) FROM
    susc_results AS s,
    rx_conv_plasma AS rxcp,
    virus_strains AS vs,
    (SELECT strain_name, COUNT(*) AS num_muts FROM strain_mutations GROUP BY strain_name) AS sm
    WHERE rxcp.ref_name = s.ref_name AND rxcp.rx_name = s.rx_name
    AND vs.strain_name = s.strain_name
    AND sm.num_muts > 1 AND sm.strain_name = s.strain_name
    AND vs.site_directed IS TRUE
-- #susc results for strains(not SDM)
SELECT COUNT(*) FROM
    susc_results AS s,
    rx_conv_plasma AS rxcp,
    virus_strains AS vs,
    (SELECT strain_name, COUNT(*) AS num_muts FROM strain_mutations GROUP BY strain_name) AS sm
    WHERE rxcp.ref_name = s.ref_name AND rxcp.rx_name = s.rx_name
    AND vs.strain_name = s.strain_name
    AND sm.strain_name = s.strain_name
    AND vs.site_directed IS FALSE
-- #IP susc results
SELECT COUNT(*) FROM susc_results AS s, rx_immu_plasma AS rxip
    WHERE rxip.ref_name = s.ref_name AND rxip.rx_name = s.rx_name;
-- #IP susc results for one mutation
SELECT COUNT(*) FROM
    susc_results AS s,
    rx_immu_plasma AS rxip,
    virus_strains AS vs,
    (SELECT strain_name, COUNT(*) AS num_muts FROM strain_mutations GROUP BY strain_name) AS sm
    WHERE rxip.ref_name = s.ref_name AND rxip.rx_name = s.rx_name
    AND vs.strain_name = s.strain_name
    AND sm.num_muts = 1 AND sm.strain_name = s.strain_name
    AND vs.site_directed IS TRUE;
-- #IP susc results for mutation combinations (not strain)
SELECT COUNT(*) FROM
    susc_results AS s,
    rx_immu_plasma AS rxip,
    virus_strains AS vs,
    (SELECT strain_name, COUNT(*) AS num_muts FROM strain_mutations GROUP BY strain_name) AS sm
    WHERE rxip.ref_name = s.ref_name AND rxip.rx_name = s.rx_name
    AND vs.strain_name = s.strain_name
    AND sm.num_muts > 1 AND sm.strain_name = s.strain_name
    AND vs.site_directed IS TRUE;
-- #IP susc results for strains (not SDM)
SELECT COUNT(*) FROM
    susc_results AS s,
    rx_immu_plasma AS rxip,
    virus_strains AS vs,
    (SELECT strain_name, COUNT(*) AS num_muts FROM strain_mutations GROUP BY strain_name) AS sm
    WHERE rxip.ref_name = s.ref_name AND rxip.rx_name = s.rx_name
    AND vs.strain_name = s.strain_name
    AND sm.strain_name = s.strain_name
    AND vs.site_directed IS FALSE;
