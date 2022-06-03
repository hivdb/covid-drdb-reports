from sql import row2dict
from preset import dump_csv
from operator import itemgetter
import math
from preset import round_number
import matplotlib.pyplot as plt
import numpy as np
from statistics import stdev, median, mean, quantiles
from mab.preset import MAIN_MAB
from preset import group_records_by
import copy
from scipy import stats
from matplotlib.ticker import MaxNLocator


log2 = lambda x: math.log(x, 2)


log10 = lambda x: math.log(x, 10)


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


def gen_omicron_mab_titer_fold(
        conn, sql, variant,
        raw_save_path,
        stat_data_path,
        figure_data_path,
        iqr_save_path):
    cursor = conn.cursor()

    cursor.execute(sql)

    records = row2dict(cursor.fetchall())
    records = adjust_potency_cmp(records)

    records = filter_and_group_records(records)

    dump_csv(raw_save_path, records)

    records = adjust_titer_and_fold_10000(records)
    records = [i for i in records if i['as_wildtype'] == 1]
    records.sort(key=itemgetter('ref_name'))

    save_results = copy.deepcopy(records)
    save_results = [
        i for i in save_results
        if not (
            i['ref_name'].startswith('FDA') and
            i['fold'] == 1 and i['ab_name'] == 'Sotrovimab')
    ]
    save_results = mark_outlier(
        save_results, 'control_ic50', 'wt_outlier',
        'wt_log_median', 'wt_log_mad')
    save_results = mark_outlier(
        save_results, 'test_ic50', 'omicron_outlier',
        'omicron_log_median', 'omicron_log_mad')
    save_results = mark_outlier(
        save_results, 'fold', 'fold_outlier',
        'fold_log_median', 'fold_log_mad')

    draw_mad_outliers(
        save_results, 'control_ic50',
        'Fold changes in IC50 relative to the normalized median IC50 of wildtype variants',
        'mAb',
        stat_data_path.parent / (
            stat_data_path.stem + '_wt' + '.png')
    )
    draw_mad_outliers(
        save_results, 'test_ic50',
        'fold change to median $IC_{50}$ (Omicron)', 'mAb',
        stat_data_path.parent / (
            stat_data_path.stem + '_omicron' + '.png')
    )
    draw_mad_outliers(
        save_results, 'fold',
        'Fold reductions in susceptibility relative to the normalized median fold reduction',
        'mAb',
        stat_data_path.parent / (
            stat_data_path.stem + '_fold' + '.png')
    )

    save_results = calc_pearsonr(save_results)

    save_results = calc_wildtype_ic50_fold(save_results)
    dump_csv(stat_data_path, save_results)

    calc_median_fold_iqr(save_results, iqr_save_path)

    dump_supplementary_table(
        copy.deepcopy(save_results), stat_data_path, variant)

    process_statistics(
        save_results, stat_data_path.parent / (
            stat_data_path.stem + '_wt' + '.csv'),
        'control_ic50',
        'ab_name', 'wt_outlier', 'wt_log_median', 'wt_log_mad')

    process_statistics(
        save_results, stat_data_path.parent / (
            stat_data_path.stem + '_omicron' + '.csv'),
        'test_ic50',
        'ab_name', 'omicron_outlier', 'omicron_log_median', 'omicron_log_mad')

    process_statistics(
        save_results, stat_data_path.parent / (
            stat_data_path.stem + '_fold' + '.csv'),
        'fold',
        'ab_name', 'fold_outlier', 'fold_log_median', 'fold_log_mad')

    # Dump figure data
    save_results = adjust_titer_and_fold_1(copy.deepcopy(records))
    save_results = [i for i in save_results if i['test_ic50']]
    dump_csv(figure_data_path, save_results)


def adjust_potency_cmp(records):
    save_results = []
    for rec in records:
        if rec['test_ic50']:

            rec['test_ic50_cmp'] = '='
            rec['control_ic50_cmp'] = '='
            if rec['test_ic50'] >= rec['test_upper_ic50']:
                rec['test_ic50_cmp'] = '>'

        del rec['test_upper_ic50']

        save_results.append(rec)

    return save_results


def skip_rec(rec):
    # Ignore because VanBalrgan tested both parent and commercial drugs
    if rec['ref_name'].startswith('VanBlargan22'):
        if rec['rx_name'] == 'AZD1061':
            return True

        if rec['rx_name'] == 'AZD7442':
            return True

        if rec['rx_name'] == 'AZD8895':
            return True

    if rec['ref_name'] == 'AstraZenica21c (FDA)':
        if rec['assay_name'] == 'Pseudovirus (FDA)':
            return True

    # if (rec['ref_name'] == 'Westendorf21'):
    #     if rec['section'] == 'Table 3D' and rec['rx_name'].endswith('_2'):
    #         return True
    #     if rec['section'] == 'Table 3B':
    #         return True

    if rec['ref_name'] == 'Cameroni21':
        if rec['assay_name'] == 'Virus isolate':
            return True

    # if (rec['ref_name'] == 'Boschi22'):
    #     return True
    # if (rec['ref_name'] == 'Gruell21'):
    #     return True

    return False


def filter_and_group_records(records):

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

            if len(ab_name_rec_list) == 1:
                results.append(ab_name_rec_list[0])
                continue

            ab_name_rec_list.sort(key=itemgetter('assay_name', 'rx_name'))
            for idx, rec in enumerate(ab_name_rec_list):
                rec['ref_name'] = '{}-{}'.format(rec['ref_name'], idx + 1)
                results.append(rec)

    return results


def adjust_titer_and_fold_10000(records):

    results = []

    for rec in records:
        if not rec['test_ic50']:
            rec['mAb'] = MAIN_MAB[rec['ab_name']]
            results.append(rec)
            continue

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

        rec['fold'] = rec['test_ic50'] / rec['control_ic50']

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
        results.append(rec)

    return results


def adjust_titer_and_fold_1(records):

    results = []

    for rec in records:
        if not rec['test_ic50']:
            rec['mAb'] = MAIN_MAB[rec['ab_name']]
            results.append(rec)
            continue

        if rec['control_ic50_cmp'] == '=' and rec['control_ic50'] < 1:
            rec['control_ic50_cmp'] == '<'
            rec['control_ic50'] = 1
        if rec['test_ic50_cmp'] == '=' and rec['test_ic50'] < 1:
            rec['test_ic50_cmp'] == '<'
            rec['test_ic50'] = 1

        rec['fold'] = rec['test_ic50'] / rec['control_ic50']

        rec['mAb'] = MAIN_MAB[rec['ab_name']]
        results.append(rec)

    return results


def short_mab_name(ab_name):
    if 'mab' in ab_name:
        return ab_name.upper()[:3]
    else:
        return ab_name


def mark_outlier(
        table, calc_column, mark_column, med_column, mad_column,
        group_by='mAb'):
    sigma = 2

    result_table = []
    for mab, mab_rec_list in group_records_by(table, group_by).items():
        values = [r[calc_column] for r in mab_rec_list]
        values = [i for i in values if i]
        median_value = median(values)
        if len(values) >= 4:
            iqr = quantiles(values, method='inclusive')
        else:
            iqr = []
        log_values = [
            log10(r[calc_column])
            for r in mab_rec_list
            if r[calc_column]
        ]
        log_median_value = median(log_values)
        log_mad = calc_mad(log_values)
        for rec in mab_rec_list:
            rec['median_{}'.format(calc_column)] = median_value
            rec['min_{}'.format(calc_column)] = min(values)
            rec['max_{}'.format(calc_column)] = max(values)
            rec['iqr_25_{}'.format(calc_column)] = iqr[0] if iqr else ''
            rec['iqr_75_{}'.format(calc_column)] = iqr[-1] if iqr else ''
            rec[med_column] = log_median_value
            rec[mad_column] = log_mad
            rec[f'log_{calc_column}'] = (
                log10(rec[calc_column])
                if rec[calc_column] else ''
            )
            rec[f'abs_dev_log_{calc_column}'] = (
                abs(log10(rec[calc_column]) - log_median_value)
                if rec[calc_column] else ''
            )
            rec[f'log2_fold_to_median_{calc_column}'] = (
                log2(float(rec[calc_column]) / median_value)
                if rec[calc_column] else ''
            )
            if not rec[calc_column]:
                rec[mark_column] = 0
            # elif is_mad_outlier(
            #         rec[calc_column], log_median_value, sigma, log_mad):
            #     rec[mark_column] = 1
            elif is_median_outlier(float(rec[calc_column]), median_value, 4):
                rec[mark_column] = 1
            else:
                rec[mark_column] = 0
            result_table.append(rec)

    return result_table


def is_median_outlier(value, median_value, fold):
    if abs(log2(value) - log2(median_value)) > log2(fold):
        return True
    else:
        return False


def is_mad_outlier(value, log_median, sigma, log_mad):
    if (log10(value) <
            log_median - sigma * log_mad):
        return True
    elif (log10(value) >
            log_median + sigma * log_mad):
        return True
    return False


def calc_pearsonr(table):

    result_table = []
    for mab, mab_rec_list in group_records_by(table, 'mAb').items():
        wv_mab_rec_list = [
            r for r in mab_rec_list
            if not r['wt_outlier'] and not r['omicron_outlier'] and not r[
                'fold_outlier']
        ]
        control_ic50_list = [
            i['control_ic50']
            for i in wv_mab_rec_list
            if i['control_ic50']
        ]
        test_ic50_list = [
            i['test_ic50']
            for i in wv_mab_rec_list
            if i['test_ic50']
        ]
        if len(control_ic50_list) < 2:
            r, p = '', ''
        else:
            r, p = stats.pearsonr(control_ic50_list, test_ic50_list)
        for rec in mab_rec_list:
            rec['r-value'] = r
            rec['p-value'] = p
            result_table.append(rec)

    return result_table


def calc_wildtype_ic50_fold(table):

    result = []
    wildtype_ic50_fold_list = []
    for mab, mab_rec_list in group_records_by(table, 'mAb').items():
        rec = mab_rec_list[0]
        wildtype_ic50_fold = rec['max_control_ic50'] / rec['min_control_ic50']
        wildtype_ic50_fold_list.append(wildtype_ic50_fold)

        for rec in mab_rec_list:
            rec['wildtype_ic50_fold'] = wildtype_ic50_fold
            result.append(rec)

    # p25, p75 = np.percentile(wildtype_ic50_fold_list, [25, 75])

    # for rec in result:
    #     rec['median_fold_fold'] = median(wildtype_ic50_fold_list)
    #     rec['iqr25_fold_fold'] = p25
    #     rec['iqr75_fold_fold'] = p75

    return result


def calc_median_fold_iqr(table, file_name):
    ab_name_list = sorted(list(set([
        r['mAb']
        for r in table
    ])))

    results = []
    for ab_name in ab_name_list:
        ab_rec_list = [r for r in table if r['mAb'] == ab_name]
        fold_list = sorted([r['fold'] for r in ab_rec_list])
        ref_list = set([r['_ref_name'] for r in ab_rec_list])
        median_fold = median(fold_list)
        if len(fold_list) < 4:
            iqr_fold = []
        else:
            iqr_fold = np.percentile(fold_list, [25, 75])
        rec = {
            'ab_name': ab_name,
            'authorized': '',
            'median_fold': round_number(median_fold),
            'iqr_25': round_number(iqr_fold[0]) if len(iqr_fold) else '',
            'iqr_75': round_number(iqr_fold[-1]) if len(iqr_fold) else '',
            'min': round_number(fold_list[0]),
            'max': round_number(fold_list[-1]),
            'fold range': round_number(fold_list[-1] / fold_list[0]),
            '#ref': len(ref_list),
            '#results': len(fold_list),
            'ref_names': ', '.join(list(sorted(ref_list))),
            'r-value': ab_rec_list[-1]['r-value'],
            'p-value': ab_rec_list[-1]['p-value'],
        }
        results.append(rec)

    dump_csv(file_name, results)


def calc_mad(records):
    """ calc median absolute deviation"""
    median_value = median(records)
    mad = median([
        abs(i - median_value)
        for i in records
    ])
    return mad


def dump_supplementary_table(table, stat_data_path, variant):
    [
        i.update({
            'fold': round_number(i['fold'])
        })
        for i in table
    ]
    [
        i.update({
            '_ref_name': i['ref_name'],
            'control_ic50': (
                "{}{}".format(
                    i.get('control_ic50_cmp'),
                    i.get('control_ic50')
                )
                if i.get('control_ic50_cmp') and i['control_ic50_cmp'] != '='
                else i.get('control_ic50', '')
            ),
            'fold': (
                i['fold']
                if i['fold_cmp'] == '='
                else f"{i['fold_cmp']}{i['fold']}"
            ),
            'omicron_ic50': (
                "{}{}".format(
                    i.get('test_ic50_cmp'),
                    i.get('test_ic50')
                )
                if i.get('test_ic50_cmp') and i['test_ic50_cmp'] != '='
                else i.get('test_ic50', '')
            ),
        })
        for i in table
    ]

    [
        i.update({
            'control_ic50': (
                f"{i['control_ic50']}*"
                if i['control_ic50'] and i['wt_outlier']
                else i['control_ic50']
            ),
            'fold': (
                f"{i['fold']}*"
                if i['fold_outlier']
                else i['fold']
            ),
            'omicron_ic50': (
                f"{i['omicron_ic50']}*"
                if i['omicron_ic50'] and i['omicron_outlier']
                else i['omicron_ic50']
            ),
        })
        for i in table
    ]

    [
        i.update({'variants': variant})
        for i in table
    ]

    headers = [
        'ref_name',
        'section',
        'ab_name',
        'variants',
        'control_var_name',
        'control_ic50',
        'omicron_ic50',
        'fold',
    ]

    dump_csv(
        stat_data_path.parent / (stat_data_path.stem + '_supply' + '.csv'),
        table, headers=headers)


def process_statistics(
        table, file_path,
        calc_column,
        column_name,
        outlier_col,
        log_median_col,
        log_mad_col):

    table = [
        i
        for i in table
        if i[calc_column]
    ]

    results = []

    for mab, mab_rec_list in group_records_by(table, column_name).items():
        num_ref = len(set([
            i['_ref_name']
            for i in mab_rec_list
        ]))
        ic_50_list = sorted([
            i[calc_column]
            for i in mab_rec_list
        ])
        iqr25, iqr75 = np.percentile(ic_50_list, [25, 75])

        low = ic_50_list[0]
        high = ic_50_list[-1]
        fold = round_number(high / low)

        # ic_50_list_no_outlier = sorted([
        #     i[calc_column]
        #     for i in mab_rec_list
        #     if not i[outlier_col]
        # ])
        # fold_no_outlier = round_number(
        #     ic_50_list_no_outlier[-1] /
        #     ic_50_list_no_outlier[0]
        # )

        med = median(ic_50_list)

        low_outliers = [
            i
            for i in mab_rec_list
            if i[outlier_col]
            and i[calc_column] <= med
            and is_median_outlier(i[calc_column], med, 4)
        ]

        within_two_fold = [
            i
            for i in mab_rec_list
            if not is_median_outlier(i[calc_column], med, 2)
        ]
        within_four_fold = [
            i
            for i in mab_rec_list
            if not is_median_outlier(i[calc_column], med, 4)
        ]

        # low_outliers = [
        #     i
        #     for i in mab_rec_list
        #     if i[outlier_col]
        #     and i[calc_column] <= med
        # ]
        low_outliers.sort(key=itemgetter(calc_column))
        low_outliers_fmt = [
            f"{i['_ref_name']} ({round_number(i[calc_column])})"
            for i in low_outliers
        ]

        high_outliers = [
            i
            for i in mab_rec_list
            if i[outlier_col]
            and i[calc_column] >= med
            and is_median_outlier(i[calc_column], med, 4)
        ]
        # high_outliers = [
        #     i
        #     for i in mab_rec_list
        #     if i[outlier_col]
        #     and i[calc_column] >= med
        # ]
        high_outliers.sort(key=itemgetter(calc_column))
        high_outliers_fmt = [
            f"{i['_ref_name']} ({round_number(i[calc_column])})"
            for i in high_outliers
        ]

        # low_outliers_ref = sorted(set(
        #         [
        #             i['_ref_name']
        #             for i in mab_rec_list
        #             if i[outlier_col]
        #             and i[calc_column] <= med
        #         ]))

        # low_outliers_value = sorted(
        #         [
        #             i[calc_column]
        #             for i in mab_rec_list
        #             if i[outlier_col]
        #             and i[calc_column] <= med
        #         ])

        # high_outliers_ref = sorted(set(
        #         [
        #             i['_ref_name']
        #             for i in mab_rec_list
        #             if i[outlier_col]
        #             and i[calc_column] >= med
        #         ]))

        # high_outliers_value = sorted(
        #         [
        #             i[calc_column]
        #             for i in mab_rec_list
        #             if i[outlier_col]
        #             and i[calc_column] >= med
        #         ])

        # log_median = mab_rec_list[0][log_median_col]
        # log_mad = mab_rec_list[0][log_mad_col]
        rec = {
            'mab': mab,
            '#ref': num_ref,
            '#result': len(mab_rec_list),
            'median': round_number(med),
            'low': round_number(low),
            'high': round_number(high),
            'high/low': fold,
            # 'fold_no_outlier': fold_no_outlier,
            'iqr25': round_number(iqr25),
            'iqr75': round_number(iqr75),
            # '#low_out': len(low_outliers_value),
            # 'low_out_ref': ','.join(low_outliers_ref),
            # 'low_out_value': ';'.join(
            # [str(round_number(i)) for i in low_outliers_value]),
            '# 2-fold': len(within_two_fold),
            '# 4-fold': len(within_four_fold),
            '#low_out': len(low_outliers_fmt),
            'low_out': ', '.join(low_outliers_fmt),
            # '#high_out': len(high_outliers_value),
            # 'high_out_ref': ','.join(high_outliers_ref),
            # 'high_out_value': ';'.join(
            # [str(round_number(i)) for i in high_outliers_value]),
            '#high_out': len(high_outliers_fmt),
            'high_out': ', '.join(high_outliers_fmt),
            # 'log_median': log_median,
            # 'log_mad': log_mad,
        }
        results.append(rec)
    dump_csv(file_path, results)


MAB_ORDER = [
        'BAM', 'ETE', 'BAM/ETE',
        'CAS', 'IMD', 'CAS/IMD',
        'SOT',
        'CIL', 'TIX', 'CIL/TIX',
        'BEB', 'ADI', 'REG',
        'AMU', 'ROM', 'AMU/ROM',
        'C135', 'C144', 'C135/C144',
]


def draw_mad_outliers(table, column, x_label, y_label, file_path):

    mab_values, mad_values = _draw_mad_outliers_mab_value(table, column)

    x_list = []
    authorized_mab_x_list = []
    authorized_mab_no_out = []
    authorized_mab_out = []
    y_list = []
    color_list = []
    ticks = []
    labels = []

    _draw_mad_outliers_prep_draw_material(
        mab_values,
        mad_values,
        x_list,
        authorized_mab_x_list,
        authorized_mab_no_out,
        authorized_mab_out,
        y_list,
        color_list,
        ticks,
        labels,
        )

    _draw_mad_outliers_draw_scatter(
        mad_values,
        x_list,
        y_list,
        color_list,
        ticks,
        labels,
        x_label,
        y_label,
        file_path
        )

    _draw_mad_outliers_draw_bin(
        authorized_mab_x_list,
        authorized_mab_no_out,
        authorized_mab_out,
        x_label,
        file_path
        )


def draw_mad_outliers_bin(table1, table2, column, x_label, y_label, file_path):

    mab_values1, mad_values1 = _draw_mad_outliers_mab_value(table1, column)
    mab_values2, mad_values2 = _draw_mad_outliers_mab_value(table2, column)

    x_list = []
    authorized_mab_x_list = []
    authorized_mab_no_out = []
    authorized_mab_out = []
    y_list = []
    color_list = []
    ticks = []
    labels = []

    _draw_mad_outliers_prep_draw_material(
        mab_values1,
        mad_values1,
        x_list,
        authorized_mab_x_list,
        authorized_mab_no_out,
        authorized_mab_out,
        y_list,
        color_list,
        ticks,
        labels,
        )
    _draw_mad_outliers_prep_draw_material(
        mab_values2,
        mad_values2,
        x_list,
        authorized_mab_x_list,
        authorized_mab_no_out,
        authorized_mab_out,
        y_list,
        color_list,
        ticks,
        labels,
        )

    _draw_mad_outliers_draw_bin(
        authorized_mab_x_list,
        authorized_mab_no_out,
        authorized_mab_out,
        x_label,
        file_path
        )


def _draw_mad_outliers_mab_value(table, column):
    mab_values = {}
    mad_values = {}
    for mab, mab_list in group_records_by(table, 'mAb').items():
        values = [
            i[column]
            for i in mab_list
            if i[column]
        ]
        log_values = [
            log2(i)
            for i in values
        ]
        log_median = median(log_values)
        abs_log_values = [
            i - log_median
            for i in log_values
        ]

        mad = calc_mad(log_values)
        low_mad = - 2 * mad
        high_mad = 2 * mad
        mad_values[mab] = (low_mad, high_mad)

        mab_values[mab] = abs_log_values

    return mab_values, mad_values


def _draw_mad_outliers_prep_draw_material(
        mab_values,
        mad_values,
        x_list,
        authorized_mab_x_list,
        authorized_mab_no_out,
        authorized_mab_out,
        y_list,
        color_list,
        ticks,
        labels,
        ):
    for odr, mab in enumerate(MAB_ORDER):
        x_values = mab_values.get(mab, [])
        if not x_values:
            continue

        low_mad, high_mad = mad_values[mab]

        x_list.extend(x_values)
        x_values_out = [
            i
            for i in x_values
            if i < low_mad or i > high_mad
        ]
        x_values_no_out = [
            i
            for i in x_values
            if i not in x_values_out
        ]
        if odr < 11:
            authorized_mab_x_list.extend(x_values)
            authorized_mab_no_out.extend(x_values_no_out)
            authorized_mab_out.extend(x_values_out)

        y_list.extend([odr] * len(x_values))
        ticks.append(odr)
        labels.append(mab)

        color_list.extend([
            'orangered' if i < low_mad or i > high_mad else 'deepskyblue'
            for i in x_values
        ])


def _draw_mad_outliers_draw_scatter(
        mad_values,
        x_list,
        y_list,
        color_list,
        ticks,
        labels,
        x_label,
        y_label,
        file_path
        ):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(x_list, y_list, c=color_list)

    for mab, (low, high) in mad_values.items():
        y = MAB_ORDER.index(mab)
        ymin = y - 0.2
        ymax = y + 0.2
        ax.vlines(x=low, ymin=ymin, ymax=ymax, color='black')
        ax.vlines(x=high, ymin=ymin, ymax=ymax, color='black')
        ax.hlines(y=y, xmin=low, xmax=high, color='black')

    max_value = math.ceil(max(x_list))
    min_value = math.floor(min(x_list))
    bins = list(np.arange(min_value, max_value + 1, 1))

    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)

    ax.set_xticks(bins)
    ax.set_xticklabels([
        rf'$2^{{{i}}}$'
        for i in bins
    ])

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    plt.savefig(
        str(file_path),
        dpi=300,
        format='png')


def _draw_mad_outliers_draw_bin(
        authorized_mab_x_list,
        authorized_mab_no_out,
        authorized_mab_out,
        x_label,
        file_path
        ):
    fig, ax = plt.subplots(figsize=(10, 10))
    max_value = math.ceil(max(authorized_mab_x_list))
    min_value = math.floor(min(authorized_mab_x_list))
    bins = list(np.arange(min_value, max_value + 1, 0.5))

    # ax.hist([
    #     authorized_mab_no_out,
    #     authorized_mab_out,
    # ], bins=bins, stacked=True)

    ax.hist(authorized_mab_x_list, bins=bins)

    ax.set_ylabel('# Results')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    bins = list(np.arange(min_value, max_value + 1, 1))
    ax.set_xticks(bins)
    ax.set_xticklabels([
        rf'$2^{{{i}}}$'
        for i in bins
    ])
    ax.set_xlabel(x_label)

    ax.axvline(x=-2, color='black', linestyle='--')
    ax.axvline(x=2, color='black', linestyle='--')

    plt.savefig(
        str(file_path.parent / (file_path.stem + '_approved_hist.png')),
        dpi=300,
        format='png'
    )
