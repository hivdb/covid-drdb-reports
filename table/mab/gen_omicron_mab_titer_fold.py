from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import row2dict
import matplotlib.pyplot as plt
import numpy as np


SUMMARY_SQL = """
SELECT DISTINCT
    s.ref_name,
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
    s.ref_name = control_pot.ref_name
    AND
    s.rx_name = control_pot.rx_name
    AND
    s.control_iso_name = control_pot.iso_name
    AND
    s.assay_name = control_pot.assay_name

    AND
    s.ref_name = test_pot.ref_name
    AND
    s.rx_name = test_pot.rx_name
    AND
    s.iso_name = test_pot.iso_name
    AND
    s.assay_name = test_pot.assay_name

    AND
    s.control_iso_name = control_iso.iso_name
    AND
    control_iso.var_name = variants.var_name

    AND
    s.iso_name = test_iso.iso_name
    AND
    test_iso.var_name = 'Omicron'
    AND
    test_iso.pattern NOT LIKE '%R346K%'

    AND
    rx.availability IS NOT NULL
;
"""


# D3 category20
CATEGORY_COLOR = [
    '#1f77b4',
    '#aec7e8',
    '#ff7f0e',
    '#ffbb78',
    '#2ca02c',
    '#98df8a',
    '#d62728',
    '#ff9896',
    '#9467bd',
    '#c5b0d5',
    '#8c564b',
    '#c49c94',
    '#e377c2',
    '#f7b6d2',
    '#7f7f7f',
    '#c7c7c7',
    '#bcbd22',
    '#dbdb8d',
    '#17becf',
    '#9edae5',
]


def gen_omicron_mab_titer_fold(
        conn,
        csv_save_path=DATA_FILE_PATH / 'mab' / 'omicron_mab_titer_fold.csv'):
    cursor = conn.cursor()

    sql = SUMMARY_SQL

    cursor.execute(sql)

    results = row2dict(cursor.fetchall())

    for rec in results:
        rec['test_ic50_cmp'] = '='
        rec['control_ic50_cmp'] = '='
        if rec['test_ic50'] >= rec['test_upper_ic50']:
            rec['test_ic50_cmp'] = '>'
        del rec['test_upper_ic50']

        if 'Monogram' in rec['assay_name']:
            rec['ref_name'] = '{} ({})'.format(rec['ref_name'], 'Monogram')
        if 'FDA' in rec['assay_name']:
            rec['ref_name'] = '{} ({})'.format(rec['ref_name'], 'FDA')

        if '_' in rec['rx_name']:
            index = rec['rx_name'].split('_')[-1]
            if not index.isdigit():
                continue
            rec['ref_name'] = '{}-{}'.format(rec['ref_name'], index)

    dump_csv(csv_save_path, results)

    results = adjust_titer_and_fold(results)

    draw_figure(
        results,
        figure_save_path=DATA_FILE_PATH / 'mab' / 'omicron_mab.svg')


def adjust_titer_and_fold(records):

    results = []

    for rec in records:
        if rec['control_ic50_cmp'] == '>' and rec['control_ic50'] >= 1000:
            rec['control_ic50'] = 10000
        if rec['test_ic50_cmp'] == '>' and rec['test_ic50'] >= 1000:
            rec['test_ic50'] = 10000
        rec['fold'] = rec['test_ic50'] / rec['control_ic50']
        results.append(rec)

        if rec['control_ic50'] < 1:
            rec['control_ic50'] = 1
            rec['control_ic50_cmp'] = '='
        if rec['test_ic50'] < 1:
            rec['test_ic50'] = 1
            rec['test_ic50_cmp'] = '='
        if rec['fold'] < 1:
            rec['fold'] = 1
            rec['fold_cmp'] = '<'

    return results


MAB_LIST = [
    'Bamlanivimab',
    'Etesevimab',
    'Bamlanivimab/Etesevimab',
    'Casirivimab',
    'Imdevimab',
    'Casirivimab/Imdevimab',
    'Cilgavimab',
    'Tixagevimab',
    'Cilgavimab/Tixagevimab',
    'Sotrovimab',
    'Regdanvimab',
    'Adintrevimab',
    'BRII-196',
    'BRII-198',
    'BRII-196/BRII-198',
]


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

            if mab_index == (len(MAB_LIST) - 1):
                sub_axis = True
            if (mab_index + cols) >= len(MAB_LIST):
                hide_x_axis = False

            if mab_index >= len(MAB_LIST):
                draw_blank(axes[row, col])
                continue

            draw_info = get_points_and_lines(
                results, MAB_LIST[mab_index], colors_map)
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
        color = CATEGORY_COLOR[idx]
        colors_map[ref_name] = color

    return colors_map


def get_points_and_lines(records, mab, colors_map):

    x_points = []
    y_points = []
    markers = []
    lines = []
    colors = []

    for rec in records:

        if not rec['as_wildtype']:
            continue

        ab_name = rec['ab_name']
        if ab_name != mab:
            continue

        control = {'ref_name': rec['ref_name']}
        test = {'ref_name': rec['ref_name']}
        fold = {'ref_name': rec['ref_name']}
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

        for rec in [control, test, fold]:

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
        'colors': colors
    }
