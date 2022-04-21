from operator import itemgetter
from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import row2dict
from preset import group_records_by
from mab.preset import MAIN_MAB
from preset import load_csv
from statistics import median
from scipy import stats

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
    test_iso.var_name IN ('Omicron/BA.1', 'Omicron/BA.2', 'Omicron/BA.1.1')

    AND
    rx.availability IS NOT NULL
;
"""


def gen_omicron_wildtype_ic50(
        conn,
        folder=DATA_FILE_PATH / 'omicron',
        file_name='omicron_wildtype_ic50_data.csv',
        report_file_name1='omicron_wildtype_ic50_11mab_stat.csv'):
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

    mab_list = [
        'Casirivimab',
        'Etesevimab',
        'Tixagevimab',
        'Bamlanivimab',
        'Cilgavimab',
        'Imdevimab',
        'Sotrovimab',
        'Bebtelovimab',
        'Casirivimab/Imdevimab',
        'Cilgavimab/Tixagevimab',
        'Bamlanivimab/Etesevimab',
    ]

    table = [i for i in table if i['ab_name'] in mab_list]

    process_virus_assay(table, folder / report_file_name1)
    create_assay_statistics(
        table, folder / 'omicron_wildtype_ic50_assay_stat.csv',
        'assay', 'PV', 'AV')

    create_assay_statistics(
        table, folder / 'omicron_wildtype_ic50_cell_line_stat.csv',
        'cell_line', '293T', 'Vero')

    create_assay_statistics(
        table, folder / 'omicron_wildtype_ic50_ACE2_1_stat.csv',
        'receptor', 'ACE2-TMPRSS2', 'ACE2')

    create_assay_statistics(
        table, folder / 'omicron_wildtype_ic50_ACE2_2_stat.csv',
        'receptor', 'ACE2-TMPRSS2', 'TMPRSS2')


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
        selector = rec['ref_name']
        selector = selector.split('-')[0]
        if selector == 'VanBlargan22':
            section = rec['section']
            if section == 'Figure 2g':
                selector = f'{selector}-1'
            else:
                selector = f'{selector}-2'

        rec['cell_line'] = cell_line_map.get(selector, '')
        rec['receptor'] = receptor_map.get(selector, '')
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

    # if rec['ref_name'] == 'Cameroni21':
    #     if rec['assay_name'] == 'Virus isolate':
    #         return True


def process_virus_assay(table, file_path):

    pseudo_list = [
        i
        for i in table
        if i['assay'] == 'PV'
    ]

    infect_list = [
        i
        for i in table
        if i['assay'] == 'AV'
    ]

    stat_detail = []
    stat_detail.append({
        'name': '#PV',
        'value': len(pseudo_list),
    })
    stat_detail.append({
        'name': '#AV',
        'value': len(infect_list),
    })

    mab_pseudo_group = group_records_by(pseudo_list, 'ab_name')
    mab_infect_group = group_records_by(infect_list, 'ab_name')

    mab_pseudo_group_length_list = [
        len(rec_list)
        for _, rec_list in mab_pseudo_group.items()
    ]

    mab_infect_group_length_list = [
        len(rec_list)
        for _, rec_list in mab_infect_group.items()
    ]

    stat_detail.append({
        'name': 'min_#PV',
        'value': min(mab_pseudo_group_length_list),
    })

    stat_detail.append({
        'name': 'max_#PV',
        'value': max(mab_pseudo_group_length_list),
    })

    stat_detail.append({
        'name': 'min_#PV',
        'value': min(mab_pseudo_group_length_list),
    })

    stat_detail.append({
        'name': 'max_#AV',
        'value': max(mab_infect_group_length_list),
    })

    stat_detail.append({
        'name': 'min_#AV',
        'value': min(mab_infect_group_length_list),
    })

    dump_csv(file_path, stat_detail)


def create_assay_statistics(
        table, file_path,
        column,
        group_a,
        group_b,
        ):

    results = []
    for mab, mab_rec_list in group_records_by(table, 'ab_name').items():

        group_a_list = [
            i['control_ic50']
            for i in mab_rec_list
            if i[column] == group_a
        ]

        group_b_list = [
            i['control_ic50']
            for i in mab_rec_list
            if i[column] == group_b
        ]

        if group_a_list:
            group_a_median = median(group_a_list)
        else:
            group_a_median = ''
        if group_b_list:
            group_b_median = median(group_b_list)
        else:
            group_b_median = ''
        if not group_a_median or not group_b_median:
            fold = ''
            r, p = '', ''
        else:
            fold = group_b_median / group_a_median
            r, p = stats.ranksums(group_a_list, group_b_list)
        rec = {
            'mab': mab,
            f'{group_a}_median': group_a_median,
            f'{group_b}_median': group_b_median,
            'fold': fold,
            'r': r,
            'p': p,
        }
        results.append(rec)

    dump_csv(file_path, results)