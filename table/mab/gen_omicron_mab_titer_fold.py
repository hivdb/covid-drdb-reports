from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import row2dict
import matplotlib.pyplot as plt


SUMMARY_SQL = """
SELECT DISTINCT
    s.ref_name,
    rx.ab_name,
    control_iso.var_name as control_var_name,
    control_pot.potency as control_ic50,
    test_iso.var_name as test_var_name,
    test_pot.potency as test_ic50,
    test_pot.potency_upper_limit as test_upper_ic50,
    s.fold_cmp,
    s.fold,
    variants.as_wildtype
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
    rx.availability IS NOT NULL
;
"""


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

    dump_csv(csv_save_path, results)

    draw_figure(
        results,
        figure_save_path=DATA_FILE_PATH / 'mab' / 'omicron_mab.svg')


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
    fig, axes = plt.subplots(5, 3, figsize=(15, 25))

    rows = axes.shape[0]
    cols = axes.shape[1]

    for row in range(rows):
        for col in range(cols):
            draw_info = get_points_and_lines(results, MAB_LIST[row * 3 + col])
            draw_sub_figure(axes[row, col], draw_info)

    plt.savefig(str(figure_save_path), format='svg', bbox_inches='tight')


def draw_sub_figure(ax, draw_info):

    ax.set_title(draw_info['mab'])

    ax.set_xticklabels(['WT', 'Omicron', 'fold'])
    ax.set_ylim([0.05, 100000])
    ax.set_yscale('log', base=10)
    ax.set_ylabel('IC50 (ng/ml) or Fold change')

    for x_points, y_points in draw_info['lines']:
        ax.plot(x_points, y_points, 'r-')

    for idx, (xp, yp, m) in enumerate(zip(
            draw_info['x_points'],
            draw_info['y_points'],
            draw_info['markers'])):

        if m == 'filled':
            ax.scatter(
                [xp], [yp], marker='o', facecolors='b', edgecolors='b')
        else:
            ax.scatter(
                [xp], [yp], marker='o', facecolors='none', edgecolors='b')


def get_points_and_lines(records, mab):

    x_points = []
    y_points = []
    markers = []
    lines = []

    for rec in records:
        if not rec['as_wildtype']:
            continue

        ab_name = rec['ab_name']
        if ab_name != mab:
            continue

        control = {}
        test = {}
        fold = {}
        for k, v in rec.items():
            if k == 'control_ic50':
                control['type'] = 'WT'
                control['value'] = v
            if k == 'test_ic50':
                test['type'] = 'Omicron'
                test['value'] = v
            if k == 'fold':
                fold['type'] = 'fold'
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
                    if v == '=':
                        markers.append('filled')
                    else:
                        markers.append('void')

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
    }
