from operator import itemgetter
from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import row2dict
from preset import group_records_by
from mab.preset import MAIN_MAB
from preset import load_csv
from statistics import median
from .preset import process_statistics
import matplotlib.pyplot as plt
import pandas as pd
from .preset import draw_mad_outliers
import seaborn as sns
import numpy as np
from .preset import mark_outlier
from preset import round_number
# from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
from .preset import MAB_ORDER


SQL = """
SELECT DISTINCT
    s.ref_name,
    s.section,
    rx.rx_name,
    rx.ab_name,
    control_iso.var_name AS control_var_name,
    control_pot.potency AS control_ic50,
    -- control_pot.potency_upper_limit AS control_upper_ic50,
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
    AND
    variants.as_wildtype == 1
;
"""


def gen_omicron_wildtype_ic50(
        conn,
        folder=DATA_FILE_PATH / 'omicron',
        file_name='wt_ic50_data.csv',
        report_file_name1='wt_ic50_11mab_stat.csv'):
    cursor = conn.cursor()

    cursor.execute(SQL)

    table = row2dict(cursor.fetchall())

    table = [i for i in table if not skip_rec(i)]
    table = [rec for rec in table if rec['ab_name'] in MAIN_MAB.keys()]
    table = [
        i
        for i in table
        if i['ref_name'] != 'Takashita22b'
    ]
    table = [
        i
        for i in table
        if not (i['ref_name'] == 'Liu21l' and i['section'] == 'Figure 2c')
    ]

    table = [
        i
        for i in table
        if not (
            i['section'] == 'Figure 3b, 3c, Extended Data Table 2'
            and
            i['ref_name'] == 'Cameroni22'
            )
    ]

    [
        rec.update({
            'assay':
                'AV' if rec['assay_name'] == 'Virus isolate' else 'PV',
            'mAb': MAIN_MAB[rec['ab_name']]
            })
        for rec in table
    ]

    table = mark_assay_on_ref_name(table)

    table.sort(key=itemgetter('ref_name'))

    cell_line_info = load_csv(
        DATA_FILE_PATH / 'omicron' / 'omicron_cellline_info.csv')

    draw_figures_mab_and_assay(table, DATA_FILE_PATH / 'omicron')

    process_cell_line(table, cell_line_info)

    # Process outliers

    draw_mad_outliers(
        table, 'control_ic50',
        'Fold changes in $IC_{50}$ relative to the normalized median $IC_{50}$ (wildtype variants)',
        'mAb',
        DATA_FILE_PATH / 'omicron' / 'wt_ic50_mad.png')

    mark_outlier(
        table, 'control_ic50', 'wt_outlier',
        'wt_log_median', 'wt_log_mad')

    dump_csv(folder / file_name, table)

    process_ref_outlier(table, folder / 'wt_ic50_ref_mab_outlier.csv')

    process_statistics(
        table, folder / 'wt_ic50_mab_stat.csv',
        'control_ic50',
        'ab_name', 'wt_outlier',
        'wt_log_median', 'wt_log_mad')

    # remove outlier
    # table = [
    #     i
    #     for i in table
    #     if not i['wt_outlier']
    # ]
    # dump_csv(folder / (file_name + 'no_outlier.csv'), table)

    # mark_outlier(
    #     table, 'control_ic50', 'assay_outlier',
    #     'wt_log_median', 'wt_log_mad', group_by='assay')
    # process_statistics(
    #     table, folder / 'wt_ic50_assay_stat.csv',
    #     'assay', 'assay_outlier')

    # mark_outlier(
    #     table, 'control_ic50', 'cl_outlier',
    #     'wt_log_median', 'wt_log_mad', group_by='cell_line')
    # process_statistics(
    #     table, folder / 'wt_ic50_cell_line_stat.csv',
    #     'cell_line', 'cl_outlier')

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
        'Adintrevimab',
        'Regdanvimab',
    ]

    table = [i for i in table if i['ab_name'] in mab_list]

    [
        i.update({
            'cell_line': (
                'other'
                if i['cell_line'] not in ['Vero', '293T', 'Huh-7']
                else i['cell_line']
            ),
            'receptor': (
                'other'
                if not i['receptor']
                else i['receptor']
            )
        })
        for i in table
    ]

    draw_figures(table, DATA_FILE_PATH / 'omicron')

    process_virus_assay(table, folder / report_file_name1)

    # analyzing different assay
    [
        i.update({
            'w_614G': (
                'D614G'
                if i['control_var_name'].startswith('B.1')
                else 'no_D614G'
            )
        })
        for i in table
    ]

    create_assay_statistics(
        table, folder / 'wt_ic50_mab_assay.csv',
        'ab_name', 'assay')

    create_assay_statistics(
        table, folder / 'wt_ic50_mab_cell_line.csv',
        'ab_name', 'cell_line')

    _table = [
        i
        for i in table
        if i['assay'] == 'AV'
    ]
    create_assay_statistics(
        _table, folder / 'wt_ic50_mab_cell_line_AV.csv',
        'ab_name', 'cell_line')

    _table = [
        i
        for i in table
        if i['assay'] == 'PV'
    ]
    create_assay_statistics(
        _table, folder / 'wt_ic50_mab_cell_line_PV.csv',
        'ab_name', 'cell_line')

    _table = [
        i
        for i in table
        if i['assay'] == 'PV'
    ]
    create_assay_statistics(
        _table, folder / 'wt_ic50_mab_receptor_PV.csv',
        'ab_name', 'receptor')

    _table = [
        i
        for i in table
        if i['assay'] == 'AV'
    ]
    create_assay_statistics(
        _table, folder / 'wt_ic50_mab_receptor_AV.csv',
        'ab_name', 'receptor')

    create_assay_statistics(
        table, folder / 'wt_ic50_mab_receptor.csv',
        'ab_name', 'receptor')

    create_assay_statistics(
        table, folder / 'wt_ic50_mab_control.csv',
        'ab_name', 'w_614G')

    create_assay_statistics(
        table, folder / 'wt_ic50_assay_cell_line.csv',
        'assay', 'cell_line')

    create_assay_statistics(
        table, folder / 'wt_ic50_assay_receptor.csv',
        'assay', 'receptor')

    create_assay_statistics(
        table, folder / 'wt_ic50_assay_control.csv',
        'assay', 'w_614G')

    create_assay_statistics(
        table, folder / 'wt_ic50_cell_line_receptor.csv',
        'cell_line', 'receptor')

    create_assay_statistics(
        table, folder / 'wt_ic50_cell_line_control.csv',
        'cell_line', 'w_614G')

    create_assay_statistics(
        table, folder / 'wt_ic50_var_receptor_control.csv',
        'receptor', 'w_614G')

    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    # ax.scatter(xs=[
    #     i['order']
    #     for i in table
    #     if i['cell_line'] != 'other'
    # ], ys=[
    #     100 if i['cell_line'] == 'Vero' else 50
    #     for i in table
    #     if i['cell_line'] != 'other'
    # ], zs=[
    #     10 if i['w_614G'] == 'D614G' else 5
    #     for i in table
    #     if i['cell_line'] != 'other'
    # ])
    # plt.show()


def process_ref_outlier(table, file_path):

    results = []
    for ref, ref_rec_list in group_records_by(table, '_ref_name').items():
        mabs = set([
            i['ab_name']
            for i in ref_rec_list
        ])
        outlier_mabs = set([
            i['ab_name']
            for i in ref_rec_list
            if i['wt_outlier']
        ])
        results.append({
            'ref_name': ref,
            '#mabs': len(mabs),
            'mabs': ', '.join(sorted(mabs)),
            '#outlier_mabs': len(outlier_mabs),
            'outlier_mab': ','.join(sorted(outlier_mabs)),
            '%outlier_mab': round_number(len(outlier_mabs) / len(mabs))
        })

    dump_csv(file_path, results)


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


def mark_assay_on_ref_name(table):

    results = []
    for ref_name, ref_name_list in group_records_by(
            table, 'ref_name').items():

        [
            i.update({'_ref_name': ref_name})
            for i in ref_name_list
        ]
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

    return results


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
        table,
        file_path,
        group_column,
        column
        ):

    results = []

    group_names = set([i[column] for i in table])

    for mab, mab_rec_list in group_records_by(table, group_column).items():

        rec = {'mab': mab}
        for g_name in group_names:
            ic50_list = sorted([
                i['control_ic50']
                for i in mab_rec_list
                if i[column] == g_name
            ])
            if len(ic50_list) > 4:
                iqr25, iqr75 = np.percentile(ic50_list, [25, 75])
            else:
                iqr25, iqr75 = '', ''
            rec[f'{g_name}_num_result'] = len(ic50_list)
            rec[f'{g_name}_median'] = round_number(
                median(ic50_list)) if ic50_list else ''
            rec[f'{g_name}_low'] = ic50_list[0] if ic50_list else ''
            rec[f'{g_name}_high'] = ic50_list[-1] if ic50_list else ''
            rec[f'{g_name}_iqr25'] = round_number(iqr25) if iqr25 else ''
            rec[f'{g_name}_iqr75'] = round_number(iqr75) if iqr75 else ''

        results.append(rec)

    dump_csv(file_path, results)


def draw_figures(table, save_folder):

    for mab, mab_rec_list in group_records_by(table, 'ab_name').items():
        for idx, rec in enumerate(
                sorted(mab_rec_list, key=itemgetter('control_ic50'))):
            rec['order'] = idx / (len(mab_rec_list) - 1)

    fig, axes = plt.subplots(4, 1, figsize=(15, 14))

    draw_figures_assay_cell_line(axes[0], table, 'assay', 'cell_line')
    draw_figures_assay_cell_line(axes[1], table, 'assay', 'receptor')
    draw_figures_assay_cell_line(axes[2], table, 'cell_line', 'receptor')
    draw_figures_assay_cell_line(axes[3], table, 'assay', 'control_var_name')

    plt.savefig(
        str(save_folder / 'assay_figure.png'),
        format='png')


def draw_figures_assay_cell_line(ax, table, col1, col2):

    # table = [
    #     i
    #     for i in table
    #     if not i['cell_line'] == 'other'
    # ]

    [
        i.update({
            'y': '{}, {}'.format(i[col1], i[col2])
        })
        for i in table
    ]

    df = pd.DataFrame(
        data={
            'x': [
                i['order'] for i in table],
            'y': [
                i['y'] for i in table],
            'hue': [
                i[col2] for i in table]
        })
    sns.stripplot(
        x="x", y="y", data=df, ax=ax)
    ax.set_xticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


def draw_figures_mab_and_assay(table, save_folder):
    fig, ax = plt.subplots(figsize=(15, 9))

    [
        i.update({
            'xlabel': f"{i['mAb']}-{i['assay']}"
        })
        for i in table
    ]

    _table = []
    for mab, mab_list in group_records_by(table, 'mAb').items():
        pv = [
            i
            for i in mab_list
            if i['assay'] == 'PV'
        ]
        av = [
            i
            for i in mab_list
            if i['assay'] == 'AV'
        ]
        if len(pv) < 2 or len(av) < 2:
            continue
        _table.extend(mab_list)

    table = _table

    figure_mab_order = [
        i
        for i in MAB_ORDER
        if i in [
            j['mAb']
            for j in table
        ]
    ]

    _table = []
    box_plot_data = []
    median_values = []
    labels = []
    ticks = []
    p_values = []
    for o in figure_mab_order:
        _list = [
            i
            for i in table
            if i['mAb'] == o
        ]
        pv_list = [
            i
            for i in _list
            if i['assay'] == 'PV'
        ]
        av_list = [
            i
            for i in _list
            if i['assay'] == 'AV'
        ]
        [
            i.update({'x': figure_mab_order.index(o) * 2 + 1})
            for i in av_list
        ]
        [
            i.update({'x': figure_mab_order.index(o) * 2 + 1.7})
            for i in pv_list
        ]
        labels.append(av_list[0]['xlabel'])
        ticks.append(figure_mab_order.index(o) * 2 + 1)
        labels.append(pv_list[0]['xlabel'])
        ticks.append(figure_mab_order.index(o) * 2 + 1.7)

        _table.extend(av_list)
        _table.extend(pv_list)

        pv_ic50 = [
            i['control_ic50']
            for i in pv_list
        ]
        av_ic50 = [
            i['control_ic50']
            for i in av_list
        ]
        median_values.append((
            median(av_ic50), median(pv_ic50)
        ))
        box_plot_data.append(av_ic50)
        box_plot_data.append(pv_ic50)
        r, p = stats.ranksums(av_ic50, pv_ic50)
        p_values.append(round_number(p))

    fold_values = []
    for i, j in median_values:
        fold = round_number(i/j)
        fold_values.append(fold)

    ax.scatter(
        [
            i['x']
            for i in table
        ],
        [
            i['control_ic50']
            for i in table
        ],
        color=[
            'orangered' if i['assay'] == 'AV' else 'deepskyblue'
            for i in table
        ]
    )
    ax.set_yscale('log', base=10)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)

    for tick in ax.get_xticklabels():
        tick.set_rotation(90)

    # [ax.axvline(x, color = 'r', linestyle='--') for x in [1,2,3,4]]
    for i in range(len(figure_mab_order)):
        x1 = i * 2 + 1
        x2 = i * 2 + 1.7
        ax.hlines(y=7000, xmin=x1, xmax=x2, color='black')
        ax.vlines(x=x1, ymin=5000, ymax=7000, color='black')
        ax.vlines(x=x2, ymin=5000, ymax=7000, color='black')
        fold_value = fold_values[i]
        p_value = p_values[i]
        p_value_star = '*' if p_value < 0.05 else ''
        ax.text(i * 2 + 1 + 0.1, 8000, f'{fold_value}{p_value_star}')

        av_median, pv_median = median_values[i]
        ax.hlines(
            y=av_median, xmin=x1 - 0.3, xmax=x1 + 0.3,
            color='black')
        ax.hlines(
            y=pv_median, xmin=x2 - 0.3, xmax=x2 + 0.3,
            color='black')

    # ax.boxplot(box_plot_data)
    # ax.set_xticks(ticks)
    # ax.set_xticklabels(labels)

    plt.savefig(
        str(save_folder / 'mab_and_assay_group.png'),
        dpi=300,
        format='png')
