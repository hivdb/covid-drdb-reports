from collections import defaultdict
from statistics import median
from resistancy import is_partial_resistant
from resistancy import is_resistant
from resistancy import is_susc
from resistancy import round_fold


def get_fold_stat(rx_list):
    fold_list = [r for r in rx_list if r['fold']]

    num_s = sum([
        r['num_fold']
        for r in fold_list
        if is_susc(r['fold'])])
    # num_s = sum(
    #     [r['num_fold']
    #     for r in fold_list
    #     if is_susc(r['fold_cmp']r['fold'])])
    num_i = sum([
        r['num_fold']
        for r in fold_list
        if is_partial_resistant(r['fold'])])
    num_r = sum([
        r['num_fold']
        for r in fold_list
        if is_resistant(r['fold'])])

    all_fold = [[i['fold']] * i['num_fold'] for i in fold_list]
    all_fold = [i for j in all_fold for i in j if i]
    median_fold = round_fold(median(all_fold)) if all_fold else ''

    num_fold = sum([r['num_fold'] for r in rx_list] + [0])

    return num_s, num_i, num_r, median_fold, num_fold


def group_var_name(var_name):
    if var_name:
        var_name = var_name.split()[0]
        var_name = var_name.split('/')[0]
        if var_name in ['Kappa' 'Iota', 'Epsilon', 'Lambda']:
            var_name = 'KIEL'
        elif var_name not in ('Alpha', 'Beta', 'Gamma', 'Delta'):
            var_name = 'other variants'
    else:
        var_name = 'other combo mut'

    return var_name


KEY_MUTATIONS = {
    'N501Y': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'N501Y'"
        ]
    },
    '∆69/70': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name LIKE '%69-70∆'"
        ]
    },
    '∆69/70 + N501Y': {
        'iso_type': 'isolate_mutations_combo_s_mut_view',
        'filter': [
            "mut.pattern LIKE '%69-70∆+N501Y'"
        ]
    },
    '∆69/70 + N501Y + A570D': {
        'iso_type': 'isolate_mutations_combo_s_mut_view',
        'filter': [
            "mut.pattern LIKE '%69-70∆+N501Y+A570D'"
        ]
    },
    '∆69/70 + N501Y + Y453F': {
        'iso_type': 'isolate_mutations_combo_s_mut_view',
        'filter': [
            "mut.pattern LIKE '%69-70∆+N501Y+Y453F'"
        ]
    },
    '∆144': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name LIKE '%144∆'"
        ]
    },
    'E484K': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'E484K'"
        ]
    },
    'Y453F': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'Y453F'"
        ]
    },
    'L452R': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'L452R'"
        ]
    },
    'E484K + N501Y': {
        'iso_type': 'isolate_mutations_combo_s_mut_view',
        'filter': [
            "mut.pattern = 'E484K+N501Y'"
        ]
    },
    'K417N': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'K417N'"
        ]
    },
    'F490S': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'F490S'"
        ]
    },
    'S494P': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'S494P'"
        ]
    },
    'K417N + E484K + N501Y': {
        'iso_type': 'isolate_mutations_combo_s_mut_view',
        'filter': [
            "mut.pattern = 'K417N+E484K+N501Y'"
        ]
    },
    'N439K': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'N439K'"
        ]
    },
    'T478K': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'T478K'"
        ]
    },
    'T95I': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'T95I'"
        ]
    },
    'G142D': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'G142D'"
        ]
    },
    'G339D': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'G339D'"
        ]
    },
    'S373F': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'S373F'"
        ]
    },
    'N440K': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'N440K'"
        ]
    },
    'G446S': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'G446S'"
        ]
    },
    'S477N': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'S477N'"
        ]
    },
    'E484A': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'E484A'"
        ]
    },
    'Q493R': {
        'iso_type': 'isolate_mutations_single_s_mut_view',
        'filter': [
            "mut.single_mut_name = 'Q493R'"
        ]
    },
}


KEY_VARIANTS = {
    'Alpha': {
        'iso_type': 'isolate_mutations_combo_s_mut_view',
        'filter': [
            "mut.var_name = 'Alpha'"
        ]
    },
    'Beta': {
        'iso_type': 'isolate_mutations_combo_s_mut_view',
        'filter': [
             "mut.var_name = 'Beta'"
        ]
    },
    'Gamma': {
        'iso_type': 'isolate_mutations_combo_s_mut_view',
        'filter': [
            "mut.var_name = 'Gamma'"
        ]
    },
    'Epsilon': {
        'iso_type': 'isolate_mutations_combo_s_mut_view',
        'filter': [
            "mut.var_name = 'Epsilon'"
        ]
    },
    'Iota': {
        'iso_type': 'isolate_mutations_combo_s_mut_view',
        'filter': [
            "mut.var_name = 'Iota'"
        ]
    },
    'Delta': {
        'iso_type': 'isolate_mutations_combo_s_mut_view',
        'filter': [
            "mut.var_name = 'Delta'"
        ]
    },
    'Kappa': {
        'iso_type': 'isolate_mutations_combo_s_mut_view',
        'filter': [
            "mut.var_name = 'Kappa'"
        ]
    },
    'Omicron': {
        'iso_type': 'isolate_mutations_combo_s_mut_view',
        'filter': [
            "mut.var_name = 'Omicron'"
        ]
    },
}
