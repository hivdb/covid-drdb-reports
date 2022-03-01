from operator import itemgetter
import math
from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import row2dict
from preset import round_number
import matplotlib.pyplot as plt
import numpy as np
from statistics import stdev, median, mean, quantiles
from .preset import MAIN_MAB
from preset import group_records_by
import copy

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
    test_iso.var_name = 'Omicron/BA.1'

    AND
    rx.availability IS NOT NULL
;
"""


# D3 category20
# CATEGORY_COLOR = [
#     '#1f77b4',
#     '#aec7e8',
#     '#ff7f0e',
#     '#ffbb78',
#     '#2ca02c',
#     '#98df8a',
#     '#d62728',
#     '#ff9896',
#     '#9467bd',
#     '#c5b0d5',
#     '#8c564b',
#     '#c49c94',
#     '#e377c2',
#     '#f7b6d2',
#     '#7f7f7f',
#     '#c7c7c7',
#     '#bcbd22',
#     '#dbdb8d',
#     '#17becf',
#     '#9edae5',
# ]

CATEGORY_COLOR = [
    '#2E91E5',
    '#E15F99',
    '#1CA71C',
    '#FB0D0D',
    '#DA16FF',
    '#222A2A',
    '#B68100',
    '#750D86',
    "#EB663B",
    "#511CFB",
    "#00A08B",
    "#FB00D1",
    "#FC0080",
    "#B2828D",
    "#6C7C32",
    "#778AAE",
    "#862A16",
    "#A777F1",
    "#620042",
    "#1616A7",
    "#DA60CA",
    "#6C4516",
    "#0D2A63",
    "#AF0038",
]


def skip_rec(rec):
    # Ignore because VanBalrgan tested both parent and commercial drugs
    if (rec['ref_name'].startswith('VanBlargan22')
            and rec['rx_name'] == 'AZD1061'):
        return True

    if (rec['ref_name'].startswith('VanBlargan22')
            and rec['rx_name'] == 'AZD7442'):
        return True

    if (rec['ref_name'].startswith('VanBlargan22')
            and rec['rx_name'] == 'AZD8895'):
        return True

    # if (rec['ref_name'] == 'Boschi22'):
    #     return True
    # if (rec['ref_name'] == 'Gruell21'):
    #     return True

    return False


def filter_records(records):

    records = [i for i in records if not skip_rec(i)]
    records = [i for i in records if i['ab_name'] in MAIN_MAB.keys()]

    for rec in records:
        rec['_ref_name'] = copy.deepcopy(rec['ref_name'])
        if 'Monogram' in rec['assay_name']:
            rec['ref_name'] = '{} ({})'.format(rec['ref_name'], 'Monogram')
        if 'FDA' in rec['assay_name']:
            rec['ref_name'] = '{} ({})'.format(rec['ref_name'], 'FDA')

        rec['assay_group'] = (
            'AV' if rec['assay_name'] == 'Virus isolate' else 'PV'
        )

    results = []
    ref_name_groups = group_records_by(records, 'ref_name')

    for ref_name, ref_name_rec_list in ref_name_groups.items():
        ab_name_groups = group_records_by(ref_name_rec_list, 'ab_name')
        for ab_name, ab_name_rec_list in ab_name_groups.items():
            if ref_name == 'Westendorf21':
                ab_name_rec_list = [
                    i
                    for i in ab_name_rec_list
                    if not (
                        i['section'] == 'Table 3D' and
                        i['rx_name'].endswith('_2')
                        ) and not (
                            i['section'] == 'Table 3B'
                        )
                ]

            if ref_name == 'Cameroni21':
                ab_name_rec_list = [
                    i
                    for i in ab_name_rec_list
                    if i['assay_group'] != 'AV'
                ]

            if len(ab_name_rec_list) == 1:
                results.append(ab_name_rec_list[0])
                continue

            ab_name_rec_list.sort(key=itemgetter('assay_name', 'rx_name'))
            for idx, rec in enumerate(ab_name_rec_list):
                rec['ref_name'] = '{}-{}'.format(rec['ref_name'], idx + 1)
                results.append(rec)

    return results


def gen_ba_1_mab_titer_fold(
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
    save_results = mark_outlier(
        save_results, 'control_ic50', 'wt_outlier',
        'wt_log_median', 'wt_log_mad')
    save_results = mark_outlier(
        save_results, 'test_ic50', 'omicron_outlier',
        'omicron_log_median', 'omicron_log_mad')
    save_results = mark_outlier(
        save_results, 'fold', 'fold_outlier',
        'fold_log_median', 'fold_log_mad')

    dump_csv(
        (
            DATA_FILE_PATH / 'mab' /
            'omicron_ba_1_mab_titer_fold_forest_figure.csv'
        ),
        save_results)

    calc_median_iqr(
        save_results, (
            DATA_FILE_PATH / 'mab' /
            'omicron_ba_1_mab_median_iqr.csv'
        ))

    report_virus_type(save_results)
    dump_for_assay_analysis(save_results)

    get_dfplot(save_results)

    # draw_figure(
    #     save_results,
    #     figure_save_path=DATA_FILE_PATH / 'mab' / 'omicron_mab.svg')

    draw_tree_figure(save_results)


def draw_tree_figure(records):

    fig, axes = plt.subplots(
        3, 1,
        gridspec_kw={'height_ratios': [16, 15, 7]},
        figsize=(5, (16 + 15 + 7) / 3)
        )
    draw_mab_figure(records, 'BAM', axes[0], 16, 5)
    draw_mab_figure(records, 'ETE', axes[1], 15, 5)
    draw_mab_figure(records, 'BAM/ETE', axes[2], 7, 5)

    plt.savefig(
        str(DATA_FILE_PATH / 'mab' / 'omicron_tree.svg'),
        format='svg', bbox_inches='tight')


def draw_mab_figure(records, mab, ax, width, height):

    records = [i for i in records if i['mAb'] == mab]

    ax.set_xlim([0.07, 12000])
    ax.set_xscale('log', base=10)
    ax.set_xlabel('IC50 (ng/ml)')

    for rec in records:
        ax.scatter(
            [rec['control_ic50']], [rec['ref_name']], color='black')
        ax.scatter(
            [rec['test_ic50']], [rec['ref_name']], color='red')
        ax.plot(
            (rec['control_ic50'], rec['test_ic50']),
            (rec['ref_name'], rec['ref_name']),
            marker='', color='black')


def dump_for_assay_analysis(records):

    prepare_records = []

    mab_list = ['BAM', 'ETE', 'CAS', 'IMD', 'CIL', 'TIX', 'SOT']
    records = [i for i in records if i['mAb'] in mab_list]

    for rec in records:
        new_rec = {
            'ref_name': "{}_{}".format(rec['ref_name'], rec['assay_group']),
            'condition': 'WT-{}'.format(rec['mAb']),
            'ic50': rec['control_ic50'],
            'assay_group': rec['assay_group'],
            'assay': rec['assay_name'],
            'section': rec['section'],
            'rx_name': rec['rx_name'],
        }
        prepare_records.append(new_rec)

    grouped_ref_name = group_records_by(prepare_records, 'ref_name')
    single_result_ref_name = set()
    for ref_name, ref_rec_list in grouped_ref_name.items():
        if len(ref_rec_list) <= 1:
            single_result_ref_name.add(ref_name)

    prepare_records = [
        i for i in prepare_records
        if i['ref_name'] not in single_result_ref_name]

    conditions = sorted(list(set(i['condition'] for i in prepare_records)))
    color_map = {
        cond: CATEGORY_COLOR[idx]
        for idx, cond in enumerate(conditions)
    }

    ref_name_sets = set([
        (i['assay_group'], i['ref_name']) for i in prepare_records
    ])

    ref_name_order_rec = [
        {
            'assay_group': i[0],
            'ref_name': i[1],
        }
        for i in ref_name_sets
    ]

    ref_name_order_rec = sorted(
        ref_name_order_rec, key=itemgetter('assay_group', 'ref_name'))

    ref_name_to_id = {
        item['ref_name']: idx + 1
        for idx, item in enumerate(ref_name_order_rec)
    }

    save_records = []
    for rec in prepare_records:
        rec['color'] = color_map[rec['condition']][1:]
        rec['order_id'] = ref_name_to_id[rec['ref_name']]
        save_records.append(rec)

    dump_csv(
        DATA_FILE_PATH / 'mab' / 'omicron_assay_analysis.csv', save_records)


def get_dfplot(dataframe):

    colors_map = get_colors_map(dataframe)

    dfplot = []

    for mab in MAIN_MAB.keys():
        draw_info = get_points_and_lines(
            dataframe, mab, colors_map)

        for (xp, yp, m, r, c) in zip(
                draw_info['x_points'],
                draw_info['y_points'],
                draw_info['markers'],
                draw_info['ref_names'],
                draw_info['colors']):
            rec = {
                'x': xp,
                'y': yp,
                'm': m,
                'ref_name': r,
                'c': c[1:],
                'mab': mab
            }
            dfplot.append(rec)

    dump_csv(DATA_FILE_PATH / 'mab' / 'omicron_mab_titer_fold_df.csv', dfplot)

    # calc_mab_cv(dfplot)


def adjust_titer_and_fold(records):

    results = []

    for rec in records:
        if rec['control_ic50_cmp'] == '>' and rec['control_ic50'] >= 1000:
            rec['control_ic50'] = 10000
        if rec['control_ic50_cmp'] == '=' and rec['control_ic50'] > 10000:
            rec['control_ic50'] = 10000
            rec['control_ic50_cmp'] = '>'

        if rec['test_ic50_cmp'] == '>' and rec['test_ic50'] >= 1000:
            rec['test_ic50'] = 10000
        if rec['test_ic50_cmp'] == '=' and rec['test_ic50'] > 10000:
            rec['test_ic50'] = 10000
            rec['test_ic50_cmp'] = '>'

        if rec['control_ic50_cmp'] == '=' and rec['control_ic50'] < 1:
            rec['control_ic50_cmp'] == '<'
            rec['control_ic50'] = 1
        if rec['test_ic50_cmp'] == '=' and rec['test_ic50'] < 1:
            rec['test_ic50_cmp'] == '<'
            rec['test_ic50'] = 1

        rec['fold'] = rec['test_ic50'] / rec['control_ic50']
        results.append(rec)

        # if rec['control_ic50'] < 1:
        #     rec['control_ic50'] = 1
        #     rec['control_ic50_cmp'] = '='
        # if rec['test_ic50'] < 1:
        #     rec['test_ic50'] = 1
        #     rec['test_ic50_cmp'] = '='
        # if rec['fold'] < 1:
        #     rec['fold'] = 1
        #     rec['fold_cmp'] = '<'

        # if '/' not in rec['ab_name']:
        #     rec['mAb'] = short_mab_name(rec['ab_name'])
        # else:
        #     ab1, ab2 = rec['ab_name'].split('/', 1)
        #     ab1, ab2 = short_mab_name(ab1), short_mab_name(ab2)
        #     rec['mAb'] = '/'.join([ab1, ab2])
        rec['mAb'] = MAIN_MAB[rec['ab_name']]

    return results


def short_mab_name(ab_name):
    if 'mab' in ab_name:
        return ab_name.upper()[:3]
    else:
        return ab_name


def draw_figure(results, figure_save_path):
    rows = 3
    cols = 6
    fig, axes = plt.subplots(rows, cols, figsize=(18, 15))

    colors_map = get_colors_map(results)

    for row in range(rows):
        hide_x_axis = False
        hide_y_axis = False
        sub_axis = False
        if row != (rows - 1):
            hide_x_axis = True

        for col in range(cols):
            if col != 0:
                hide_y_axis = True
            if col == (cols - 1):
                sub_axis = True

            mab_index = row * cols + col

            if mab_index == (len(MAIN_MAB) - 1):
                sub_axis = True
            if (mab_index + cols) >= len(MAIN_MAB):
                hide_x_axis = False

            if mab_index >= len(MAIN_MAB):
                draw_blank(axes[row, col])
                continue

            draw_info = get_points_and_lines(
                results, list(MAIN_MAB.keys())[mab_index], colors_map)
            draw_sub_figure(
                axes[row, col], draw_info, hide_x_axis, hide_y_axis, sub_axis)

    draw_legend(axes[rows-1, cols-2], colors_map)
    fig.subplots_adjust(wspace=0, hspace=0)

    plt.savefig(str(figure_save_path), format='svg', bbox_inches='tight')


def draw_legend(ax, colors_map):

    for k, v in colors_map.items():
        ax.scatter(
                [],
                [],
                marker='o',
                c=np.array([v]),
                label=k)

    ax.legend(loc="center left", ncol=2)


def draw_blank(ax):
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.patch.set_alpha(0)


def draw_sub_figure(
        ax, draw_info, hide_x_axis=False, hide_y_axis=False, sub_axis=False):

    ax.set_title(draw_info['mab'], y=1.0, pad=-14, loc='left')

    y_lower = 0.7
    y_upper = 30000

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['WT', 'Omicron', 'Fold'])
    ax.set_ylim([y_lower, y_upper])
    ax.set_yscale('log', base=10)
    ax.set_ylabel('IC50 (ng/ml)')
    ax.yaxis.label.set_color('blue')
    ax.tick_params(axis='y', colors='blue')

    if hide_x_axis:
        ax.get_xaxis().set_visible(False)
    if hide_y_axis:
        ax.get_yaxis().set_visible(False)

    # for x_points, y_points in draw_info['lines']:
    #     ax.plot(x_points, y_points, 'k-')

    draw_points(ax, draw_info)

    if sub_axis:
        ax2 = ax.twinx()
        ax2.set_ylim([y_lower, y_upper])
        ax2.set_yscale('log', base=10)
        ax2.set_ylabel('Fold change')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')

        draw_points(ax2, draw_info)


def draw_points(ax, draw_info):
    wt_x_points = [
        i for idx, i in enumerate(draw_info['y_points']) if (idx + 1) % 3 == 1]
    omicron_x_points = [
        i for idx, i in enumerate(draw_info['y_points']) if (idx + 1) % 3 == 2]
    fold_x_points = [
        i for idx, i in enumerate(draw_info['y_points']) if (idx + 1) % 3 == 0]
    ax.boxplot(
        [wt_x_points, omicron_x_points, fold_x_points],
        positions=[0, 1, 2],
        labels=['WT', 'Omicron', 'Fold'],
        showfliers=False,
        boxprops={'color': 'gray'},
        capprops={'color': 'gray'},
        whiskerprops={'color': 'gray'})

    for idx, (xp, yp, m, c) in enumerate(zip(
            draw_info['x_points'],
            draw_info['y_points'],
            draw_info['markers'],
            draw_info['colors'])):

        ax.scatter(
                [xp], [yp], marker=m, facecolors=c, edgecolors=c)


def get_colors_map(records):
    ref_list = sorted(set([rec['ref_name'] for rec in records]))

    colors_map = {}
    for idx, ref_name in enumerate(ref_list):
        if idx >= len(CATEGORY_COLOR):
            color = '#000000'
        else:
            color = CATEGORY_COLOR[idx]
        colors_map[ref_name] = color

    return colors_map


def get_points_and_lines(records, mab, colors_map):

    x_points = []
    y_points = []
    markers = []
    lines = []
    colors = []
    ref_names = []

    for rec in records:

        ab_name = rec['ab_name']
        if ab_name != mab:
            continue

        ref_name = rec['ref_name']

        control = {'ref_name': ref_name}
        test = {'ref_name': ref_name}
        fold = {'ref_name': ref_name}
        for k, v in rec.items():
            if k == 'control_ic50':
                control['type'] = 0
                control['value'] = v
            if k == 'test_ic50':
                test['type'] = 1
                test['value'] = v
            if k == 'fold':
                fold['type'] = 2
                fold['value'] = v
            if k == 'control_ic50_cmp':
                control['cmp'] = v
            if k == 'test_ic50_cmp':
                test['cmp'] = v
            if k == 'fold_cmp':
                fold['cmp'] = v

        if not rec['as_wildtype']:
            rec_list = [test]
        else:
            rec_list = [control, test, fold]

        for rec in rec_list:

            for k, v in rec.items():
                if k == 'type':
                    x_points.append(v)
                if k == 'value':
                    y_points.append(v)
                if k == 'cmp':
                    if v == '>':
                        markers.append('^')
                    # elif v == '<':
                    #     markers.append('v')
                    else:
                        markers.append('o')

            ref_name = rec['ref_name']
            choose_color = colors_map[ref_name]
            colors.append(choose_color)
            ref_names.append(ref_name)

        lines.append((
            x_points[-3:-1],
            y_points[-3:-1]
        ))
        lines.append((
            x_points[-2:],
            y_points[-2:]
        ))

    return {
        'mab': mab,
        'lines': lines,
        'x_points': x_points,
        'y_points': y_points,
        'markers': markers,
        'colors': colors,
        'ref_names': ref_names
    }


def calc_coefficient_of_variation(dataframe, log_sample=False):
    samples = [i['y'] for i in dataframe]

    if log_sample:
        samples = [math.log(i) for i in samples]

    mean_value = mean(samples)
    std = stdev(samples)

    return std / mean_value


def calc_mab_cv(records):

    mabs = list(set([i['mab'] for i in records]))

    results = []
    for mab in mabs:
        for idx, t in zip([0, 1, 2], ['wt', 'omicron', 'fold']):
            dataframe = [
                i for i in records if i['mab'] == mab and i['x'] == idx]
            if len(dataframe) > 1:
                cv_value = calc_coefficient_of_variation(dataframe)
                cv_value = (cv_value * 100 // 1) / 100
                cv_value_log = calc_coefficient_of_variation(
                    dataframe, log_sample=True)
                cv_value_log = (cv_value_log * 100 // 1) / 100
            else:
                cv_value = 'NA'
                cv_value_log = 'NA'
            results.append({
                'mab': mab,
                'data_type': t,
                'cv': cv_value,
                'cv_log_sample': cv_value_log,
                'num_sample': len(dataframe),
                # 'samples': [math.log(i['y']) for i in dataframe]
            })

    dump_csv(DATA_FILE_PATH / 'mab' / 'omicron_mab_cv.csv', results)


def report_virus_type(records):
    report_records = []
    mab_group = group_records_by(records, 'ab_name')
    for mab, mab_rec_list in mab_group.items():
        assay_group = group_records_by(mab_rec_list, 'assay_group')
        report_rec = {
            'mab': mab
        }
        for assay, assay_list in assay_group.items():
            median_ic50 = median([i['control_ic50'] for i in assay_list])
            report_rec[assay] = median_ic50
            report_rec['num_{}'.format(assay)] = len(assay_list)
        if 'AV' in report_rec and 'PV' in report_rec and report_rec['PV'] > 0:
            report_rec['fold'] = report_rec['AV'] / report_rec['PV']
        else:
            report_rec['fold'] = ''
        report_records.append(report_rec)

    dump_csv(
        DATA_FILE_PATH / 'mab' / 'omicron_assay_virus_type_report.csv',
        report_records
    )


def calc_median_iqr(table, file_name):
    ab_name_list = sorted(list(set([
        r['mAb']
        for r in table
    ])))

    results = []
    for ab_name in ab_name_list:
        fold_list = sorted([r['fold'] for r in table if r['mAb'] == ab_name])
        ref_list = set([r['ref_name'] for r in table if r['mAb'] == ab_name])
        median_fold = median(fold_list)
        if len(fold_list) < 4:
            iqr_fold = []
        else:
            iqr_fold = quantiles(fold_list, method='inclusive')
        rec = {
            'ab_name': ab_name,
            'authorized': '',
            'median_fold': round_number(median_fold),
            'iqr_25': round_number(iqr_fold[0]) if iqr_fold else '',
            'iqr_75': round_number(iqr_fold[-1]) if iqr_fold else '',
            'min': round_number(fold_list[0]),
            'max': round_number(fold_list[-1]),
            'fold range': round_number(fold_list[-1] / fold_list[0]),
            'num_ref': len(ref_list),
            'ref_names': ', '.join(list(sorted(ref_list)))
        }
        results.append(rec)

    dump_csv(file_name, results)


def get_median_absolute_deviation(records):
    median_value = median(records)
    mad = median([
        abs(i - median_value)
        for i in records
    ])
    return mad


def mark_outlier(table, calc_column, mark_column, med_column, mad_column):
    sigma = 2

    result_table = []
    for mab, mab_rec_list in group_records_by(table, 'mAb').items():
        values = [r[calc_column] for r in mab_rec_list]
        median_value = median(values)
        if len(values) >= 4:
            iqr = quantiles(values, method='inclusive')
        else:
            iqr = []
        log_values = [
            math.log(r[calc_column], 10)
            for r in mab_rec_list
        ]
        log_median_value = median(log_values)
        log_mad = get_median_absolute_deviation(log_values)
        for rec in mab_rec_list:
            rec['median_{}'.format(calc_column)] = median_value
            rec['iqr_25_{}'.format(calc_column)] = iqr[0] if iqr else ''
            rec['iqr_75_{}'.format(calc_column)] = iqr[-1] if iqr else ''
            rec[med_column] = log_median_value
            rec[mad_column] = log_mad
            if (
                    math.log(rec[calc_column], 10) <
                    log_median_value - sigma * log_mad):
                rec[mark_column] = 1
            elif (
                    math.log(rec[calc_column], 10) >
                    log_median_value + sigma * log_mad):
                rec[mark_column] = 1
            else:
                rec[mark_column] = 0
            result_table.append(rec)

    return result_table
