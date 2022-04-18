from operator import itemgetter
from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import row2dict
from preset import group_records_by
from mab.preset import MAIN_MAB
from preset import load_csv

SQL = """
SELECT DISTINCT
    s.ref_name,
    s.section,
    rx.rx_name,
    rx.ab_name,
    control_iso.var_name AS control_var_name,
    control_pot.potency AS control_ic50,
    control_pot.potency_upper_limit AS control_upper_ic50,
    variants.as_wildtype,
    s.assay_name
FROM
    susc_results_view s,
    rx_mab_view rx,
    rx_potency control_pot,
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
    s.control_iso_name = control_iso.iso_name
    AND
    control_iso.var_name = variants.var_name

    AND
    s.iso_name = test_iso.iso_name
    AND
    test_iso.var_name IN ('Omicron/BA.1', 'Omicron/BA.2', 'Omicron/BA.1.1')

    AND
    rx.availability IS NOT NULL
;
"""


def gen_omicron_wildtype_ic50(
        conn,
        folder=DATA_FILE_PATH / 'omicron',
        file_name='omicron_wildtype_ic50_data.csv'):
    cursor = conn.cursor()

    cursor.execute(SQL)

    table = row2dict(cursor.fetchall())

    table = [i for i in table if not skip_rec(i)]
    table = [rec for rec in table if rec['ab_name'] in MAIN_MAB.keys()]

    [
        rec.update({
            'assay':
                'AV' if rec['assay_name'] == 'Virus isolate' else 'PV',
            'mAb': MAIN_MAB[rec['ab_name']]
            })
        for rec in table
    ]

    results = []
    for ref_name, ref_name_list in group_records_by(
            table, 'ref_name').items():
        pv_list = [
            i
            for i in ref_name_list
            if i['assay'] == 'PV'
        ]

        av_list = [
            i
            for i in ref_name_list
            if i['assay'] == 'AV'
        ]

        if not pv_list:
            for rec in av_list:
                results.append(rec)
        elif not av_list:
            for rec in pv_list:
                results.append(rec)
        else:
            for rec in av_list:
                rec['ref_name'] = f'{ref_name}-1'
                results.append(rec)
            for rec in pv_list:
                rec['ref_name'] = f'{ref_name}-2'
                results.append(rec)

    results.sort(key=itemgetter('ref_name'))

    cell_line_info = load_csv(
        DATA_FILE_PATH / 'omicron' / 'omicron_cellline_info.csv')

    process_cell_line(results, cell_line_info)

    dump_csv(folder / file_name, results)


def process_cell_line(table, cell_line_info):
    cell_line_map = {
        rec['Ref']: rec['cell']
        for rec in cell_line_info
    }

    receptor_map = {
        rec['Ref']: rec['receptor']
        for rec in cell_line_info
    }

    results = []
    for rec in table:
        rec['cell_line'] = cell_line_map.get(rec['ref_name'], '')
        rec['receptor'] = receptor_map.get(rec['ref_name'], '')
        results.append(rec)

    return results


def skip_rec(rec):
    # Ignore because VanBalrgan tested both parent and commercial drugs
    if rec['ref_name'].startswith('VanBlargan22'):
        if rec['rx_name'] == 'AZD1061':
            return True

        if rec['rx_name'] == 'AZD7442':
            return True

        if rec['rx_name'] == 'AZD8895':
            return True

    # if (rec['ref_name'] == 'Westendorf21'):
    #     if rec['section'] == 'Table 3D' and rec['rx_name'].endswith('_2'):
    #         return True
    #     if rec['section'] == 'Table 3B':
    #         return True

    if rec['ref_name'] == 'Cameroni21':
        if rec['assay_name'] == 'Virus isolate':
            return True
