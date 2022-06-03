from .preset import skip_rec
from preset import DATA_FILE_PATH
from preset import dump_csv
from sql import row2dict

SUMMARY_SQL = """
SELECT DISTINCT
    ba1.ref_name,
    ba1.section,
    ba1.potency_type,
    ba1.control_iso_name AS control_var_name,
    ba1.control_potency,
    rx.rx_name,
    rx.ab_name,
    ba1.fold_cmp AS ba1_fold_cmp,
    ba1.fold AS ba1_fold,
    ba1.potency AS ba1_ic50,
    ba11.fold_cmp AS ba11_fold_cmp,
    ba11.fold AS ba11_fold,
    ba11.potency AS ba11_ic50,
    ba1.assay_name
FROM
    susc_results_view ba1,
    susc_results_view ba11,
    rx_mab_view rx,
    isolate_mutations_combo_s_mut_view ba1_iso,
    isolate_mutations_combo_s_mut_view ba11_iso
WHERE
    ba1.ref_name = ba11.ref_name
    AND
    ba1.rx_name = ba11.rx_name
    AND
    ba1.control_iso_name = ba11.control_iso_name
    AND
    ba1.potency_type = ba11.potency_type
    AND
    ba1.potency_type = 'IC50'

    AND
    ba1.ref_name = rx.ref_name
    AND
    ba1.rx_name = rx.rx_name

    AND
    ba1.iso_name = ba1_iso.iso_name
    AND
    ba1_iso.var_name = 'Omicron/BA.1'

    AND
    ba11.iso_name = ba11_iso.iso_name
    AND
    ba11_iso.var_name = 'Omicron/BA.1.1'

    AND
    rx.availability IS NOT NULL
;
"""


def gen_omicron_titer_fold_ba1_compare(
        conn,
        folder=DATA_FILE_PATH / 'omicron',
        file_name='omicron_ba1_compare.csv'):
    cursor = conn.cursor()

    cursor.execute(SUMMARY_SQL)

    table = row2dict(cursor.fetchall())

    table = [i for i in table if not skip_rec(i)]

    dump_csv(folder / file_name, table)
