from preset import row2dict
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
        conn, sql,
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
        if not (i['ref_name'].startswith('FDA') and i['fold'] == 1 and i['ab_name'] == 'Sotrovimab')
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

    save_results = calc_pearsonr(save_results)

    save_results = calc_wildtype_ic50_fold(save_results)
    dump_csv(stat_data_path, save_results)

    calc_median_fold_iqr(save_results, iqr_save_path)

    headers = [
        'ref_name',
        'section',
        'rx_name',
        'ab_name',
        'control_var_name',
        'control_ic50',
        'test_var_name',
        'test_ic50',
        'fold_cmp',
        'fold',
        'as_wildtype',
        'assay_name',
        'test_ic50_cmp',
        'control_ic50_cmp',
        '_ref_name',
        'assay_group',
        'mAb',
        'median_control_ic50',
        'min_control_ic50',
        'max_control_ic50',
        'iqr_25_control_ic50',
        'iqr_75_control_ic50',
        'wt_log_median',
        'wt_log_mad',
        'wt_outlier',
    ]
    dump_csv(
        figure_data_path.parent / (figure_data_path.stem + '_0' + '.csv'),
        save_results, headers=headers)

    headers = [
        'ref_name',
        'section',
        'rx_name',
        'ab_name',
        'control_var_name',
        'control_ic50',
        'test_var_name',
        'test_ic50',
        'fold_cmp',
        'fold',
        'as_wildtype',
        'assay_name',
        'test_ic50_cmp',
        'control_ic50_cmp',
        '_ref_name',
        'assay_group',
        'mAb',
        'median_test_ic50',
        'min_test_ic50',
        'max_test_ic50',
        'iqr_25_test_ic50',
        'iqr_75_test_ic50',
        'omicron_log_median',
        'omicron_log_mad',
        'omicron_outlier',
    ]
    dump_csv(
        figure_data_path.parent / (figure_data_path.stem + '_1' + '.csv'),
        save_results, headers=headers)

    headers = [
        'ref_name',
        'section',
        'rx_name',
        'ab_name',
        'control_var_name',
        'control_ic50',
        'test_var_name',
        'test_ic50',
        'fold_cmp',
        'fold',
        'as_wildtype',
        'assay_name',
        'test_ic50_cmp',
        'control_ic50_cmp',
        '_ref_name',
        'assay_group',
        'mAb',
        'median_fold',
        'min_fold',
        'max_fold',
        'iqr_25_fold',
        'iqr_75_fold',
        'fold_log_median',
        'fold_log_mad',
        'fold_outlier',
        'r-value',
        'p-value',
    ]
    dump_csv(
        figure_data_path.parent / (figure_data_path.stem + '_2' + '.csv'),
        save_results, headers=headers)

    headers = [
        'ref_name',
        'section',
        'rx_name',
        'ab_name',
        'control_var_name',
        'control_ic50',
        'test_var_name',
        'test_ic50',
        'fold_cmp',
        'fold',
        'as_wildtype',
        'assay_name',
        'test_ic50_cmp',
        'control_ic50_cmp',
        '_ref_name',
        'assay_group',
        'mAb',
        'wildtype_ic50_fold',
    ]
    dump_csv(
        figure_data_path.parent / (figure_data_path.stem + '_3' + '.csv'),
        save_results, headers=headers)

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

    if rec['ref_name'] == 'Cao22':
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


def mark_outlier(table, calc_column, mark_column, med_column, mad_column):
    sigma = 2

    result_table = []
    for mab, mab_rec_list in group_records_by(table, 'mAb').items():
        values = [r[calc_column] for r in mab_rec_list]
        values = [i for i in values if i]
        median_value = median(values)
        if len(values) >= 4:
            iqr = quantiles(values, method='inclusive')
        else:
            iqr = []
        log_values = [
            math.log(r[calc_column], 10)
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
            if not rec[calc_column]:
                rec[mark_column] = 0
            elif (
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
        fold_list = sorted([r['fold'] for r in table if r['mAb'] == ab_name])
        ref_list = set([r['_ref_name'] for r in table if r['mAb'] == ab_name])
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
            'num_ref': len(ref_list),
            'ref_names': ', '.join(list(sorted(ref_list)))
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
