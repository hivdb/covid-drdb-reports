from collections import defaultdict
from statistics import median
from resistancy import is_partial_resistant
from resistancy import is_resistant
from resistancy import is_susc
from resistancy import round_fold
from pathlib import Path
from preset import load_yaml


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
    ret_var_name = var_name
    if ret_var_name:
        ret_var_name = ret_var_name.split()[0]
        if not ret_var_name.startswith('Omicron'):
            ret_var_name = ret_var_name.split('/')[0]
        if ret_var_name in [
                'Kappa' 'Iota',
                'Epsilon', 'Lambda', 'Eta', 'Mu']:
            ret_var_name = 'VOIs (Kappa, Iota, Epsilon, Lambda, Eta, Mu)'
        elif ret_var_name not in (
                'Alpha', 'Beta', 'Gamma', 'Delta') and (
                    not ret_var_name.startswith('Omicron')
                ):
            ret_var_name = 'other variants'
    else:
        ret_var_name = 'other combo mut'
    return ret_var_name


def load_key_mutations(mapper, key_mut_file):

    for mut, match in load_yaml(key_mut_file).items():
        if '+' in match:
            table_name = 'isolate_mutations_combo_s_mut_view'
            if '%' in match:
                filter = f"mut.pattern LIKE '{match}'"
            else:
                filter = f"mut.pattern = '{match}'"
        else:
            table_name = 'isolate_mutations_single_s_mut_view'
            if '%' in match:
                filter = f"mut.single_mut_name LIKE '{match}'"
            else:
                filter = f"mut.single_mut_name = '{match}'"

        mapper[mut] = {
            'iso_type': table_name,
            'filter': [
                filter
            ]
        }


def load_key_variants(mapper, var_file):

    for mut, match in load_yaml(var_file).items():
        table_name = 'isolate_mutations_combo_s_mut_view'

        if type(match) == list:
            match = [
                f'"{i}"'
                for i in match
            ]
            filter = f"mut.var_name IN ({', '.join(match)})"
        elif '%' in match:
            filter = f"mut.var_name LIKE '{match}'"
        else:
            filter = f"mut.var_name = '{match}'"

        mapper[mut] = {
            'iso_type': table_name,
            'filter': [
                filter
            ]
        }


KEY_MUTATIONS = {}
load_key_mutations(
    KEY_MUTATIONS,
    Path(__file__).resolve().parent / 'key_mutations.yml')


OMICRON_MUTATIONS = {}

load_key_mutations(
    OMICRON_MUTATIONS,
    Path(__file__).resolve().parent / 'omicron_mutations.yml')


KEY_VARIANTS = {}
load_key_variants(
    KEY_VARIANTS,
    Path(__file__).resolve().parent / 'key_variants.yml')


OMICRON_VARIANTS = {}
load_key_variants(
    OMICRON_VARIANTS,
    Path(__file__).resolve().parent / 'omicron_variants.yml')
