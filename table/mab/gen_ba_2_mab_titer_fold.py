from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import row2dict
from .gen_omicron_mab_titer_fold import adjust_titer_and_fold
from .gen_omicron_mab_titer_fold import filter_records

SUMMARY_SQL = """
SELECT DISTINCT
    s.ref_name,
    s.section,
    rx.rx_name,
    rx.ab_name,
    control_iso.var_name as control_var_name,
    control_pot.potency as control_ic50,
    test_iso.var_name as test_var_name,
    test_pot.potency as test_ic50,
    test_pot.potency_upper_limit as test_upper_ic50,
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
    test_iso.var_name = 'Omicron/BA.2'

    AND
    rx.availability IS NOT NULL
;
"""


def gen_ba_2_mab_titer_fold(
        conn,
        csv_save_path=DATA_FILE_PATH / 'mab' / 'omicron_mab_titer_fold.csv'):
    cursor = conn.cursor()

    sql = SUMMARY_SQL

    cursor.execute(sql)

    records = row2dict(cursor.fetchall())

    records = filter_records(records)

    save_results = []
    for rec in records:
        rec['test_ic50_cmp'] = '='
        rec['control_ic50_cmp'] = '='
        if rec['test_ic50'] >= rec['test_upper_ic50']:
            rec['test_ic50_cmp'] = '>'
        del rec['test_upper_ic50']

        save_results.append(rec)

    dump_csv(csv_save_path, save_results)

    save_results = adjust_titer_and_fold(save_results)
    save_results = [i for i in save_results if i['as_wildtype'] == 1]

    dump_csv(
        DATA_FILE_PATH / 'mab' / 'omicron_ba_2_mab_titer_fold_forest_figure.csv',
        save_results)

