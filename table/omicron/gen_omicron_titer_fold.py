from .preset import gen_omicron_mab_titer_fold
from preset import DATA_FILE_PATH
from .preset import draw_mad_outliers_bin
from .preset import draw_mad_outliers
from preset import load_csv

SUMMARY_SQL = """
SELECT DISTINCT
    s.ref_name,
    s.section,
    rx.rx_name,
    rx.ab_name,
    control_iso.var_name AS control_var_name,
    control_pot.potency AS control_ic50,
    test_iso.var_name AS test_var_name,
    test_pot.potency AS test_ic50,
    test_pot.potency_upper_limit AS test_upper_ic50,
    s.fold_cmp,
    s.fold,
    variants.as_wildtype,
    s.assay_name
FROM
    susc_results_view s,
    rx_mab_view rx,
    rx_potency control_pot,
    rx_potency test_pot,
    isolates control_iso,
    variants,
    isolate_mutations_combo_s_mut_view test_iso
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.potency_type = 'IC50'

    AND
    s.ref_name = control_pot.ref_name
    AND
    s.rx_name = control_pot.rx_name
    AND
    s.control_iso_name = control_pot.iso_name
    AND
    s.assay_name = control_pot.assay_name
    AND
    s.potency_type = control_pot.potency_type

    AND
    s.ref_name = test_pot.ref_name
    AND
    s.rx_name = test_pot.rx_name
    AND
    s.iso_name = test_pot.iso_name
    AND
    s.assay_name = test_pot.assay_name
    AND
    s.potency_type = test_pot.potency_type

    AND
    s.control_iso_name = control_iso.iso_name
    AND
    control_iso.var_name = variants.var_name

    AND
    s.iso_name = test_iso.iso_name
    AND
    test_iso.var_name = '{variant}'

    AND
    rx.availability IS NOT NULL

UNION

SELECT DISTINCT
    s.ref_name,
    s.section,
    rx.rx_name,
    rx.ab_name,
    s.control_iso_name AS control_var_name,
    '' AS control_ic50,
    s.iso_name AS test_var_name,
    '' AS test_ic50,
    '' AS test_upper_ic50,
    s.fold_cmp,
    s.fold,
    variants.as_wildtype,
    s.assay_name
FROM
    susc_results_view s,
    rx_mab_view rx,
    isolates control_iso,
    variants,
    isolate_mutations_combo_s_mut_view test_iso
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.potency_type = 'IC50'

    AND
    s.potency IS NULL

    AND
    s.control_iso_name = control_iso.iso_name
    AND
    control_iso.var_name = variants.var_name

    AND
    s.iso_name = test_iso.iso_name
    AND
    test_iso.var_name = '{variant}'

    AND
    rx.availability IS NOT NULL
;
"""


def gen_omicron_titer_fold(conn):
    omicron_path = DATA_FILE_PATH / 'omicron'

    for variant in ['Omicron/BA.1', 'Omicron/BA.2', 'Omicron/BA.1.1']:
        sql = SUMMARY_SQL.format(variant=variant)
        variant_name = variant.split('/', 1)[-1]
        gen_omicron_mab_titer_fold(
            conn, sql,
            variant,
            omicron_path / f'{variant_name}_raw_data.csv',
            omicron_path / f'{variant_name}_stat_data.csv',
            omicron_path / f'{variant_name}_figure_data.csv',
            omicron_path / f'{variant_name}_iqr_data.csv')

    table_list = []
    table_list2 = []
    for variant in ['Omicron/BA.1', 'Omicron/BA.2']:
        sql = SUMMARY_SQL.format(variant=variant)
        variant_name = variant.split('/', 1)[-1]
        table = load_csv(omicron_path / f'{variant_name}_stat_data.csv')
        [
            i.update({'fold': float(i['fold'])})
            for i in table
        ]
        table_list.append(table)
        table_list2.extend(table)

    draw_mad_outliers_bin(
        table_list[0], table_list[1], 'fold',
        'Fold reductions in susceptibility relative to the normalized median fold reduction',
        'mAb',
        omicron_path / 'BA_1_BA_2_fold.png')

    draw_mad_outliers(
        table_list2, 'fold',
        'Fold reductions in susceptibility relative to the normalized median fold reduction',
        'mAb',
        omicron_path / 'BA_1_BA_2_fold_merge_median.png')
