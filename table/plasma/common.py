import statistics
from operator import itemgetter
from collections import defaultdict
from preset import round_number
from preset import dump_json
from resistancy import RESISTANCE_FILTER
from resistancy import is_susc
from resistancy import is_partial_resistant
from resistancy import is_resistant
from statistics import median
from .preset import PLASMA_RENAME
from .preset import PLASMA_POST_RENAME
from .preset import RENAME_CP_EXECUTOR


def gen_plasma_indiv_table(
        conn, row_filters, subrow_filters,
        sql_template, record_modifier=None):

    cursor = conn.cursor()

    records = defaultdict(dict)
    for row_name, attr_r in row_filters.items():
        for subrow_name, attr_subr in subrow_filters.items():
            for resist_name, resist_filter in RESISTANCE_FILTER.items():
                rxtype = attr_subr['rxtype']

                r_filter = attr_r.get('filter', [])
                filter = '\n    '.join(r_filter + resist_filter)

                if subrow_name.lower().startswith('cp'):
                    filter += '\n   '
                    filter += '\n   '.join(attr_subr.get('cp_filters', []))

                sql = sql_template.format(
                    rxtype=rxtype,
                    filters=filter
                )
                # print(sql)

                cursor.execute(sql)
                for row in cursor.fetchall():
                    variant_name = row_name
                    cp_name = row['rx_name']
                    reference = row['ref_name']
                    num_results = row['sample_count']
                    fold = row['fold']

                    # if cp_name in EXCLUDE_PLASMA:
                    #     continue
                    # exclude_tester = EXCLUDE_STUDIES.get(reference)
                    # if exclude_tester and exclude_tester(cp_name):
                    #     continue

                    cp_name = PLASMA_RENAME.get(cp_name, cp_name)
                    rename_executors = RENAME_CP_EXECUTOR.get(reference, [])
                    for tester, new_name in rename_executors:
                        if tester(cp_name):
                            cp_name = new_name

                    key = '{}{}{}'.format(variant_name, cp_name, reference)

                    rec = records[key]
                    rec['Variant name'] = variant_name
                    rec['Plasma'] = PLASMA_POST_RENAME.get(cp_name, cp_name)
                    rec['S'] = rec.get('S', 0)
                    rec['I'] = rec.get('I', 0)
                    rec['R'] = rec.get('R', 0)

                    if not rec.get('folds'):
                        rec['folds'] = []

                    fold_list = rec.get('folds')
                    if fold is not None:
                        fold_list.append(fold)

                    if resist_name == 'susceptible':
                        rec['S'] += num_results
                    elif resist_name == 'partial':
                        rec['I'] += num_results
                    else:
                        rec['R'] += num_results

                    rec['Reference'] = reference

    for rec in records.values():
        folds = rec.get('folds', [])
        if folds:
            rec['Median'] = str(round_number(statistics.median(folds)))
        else:
            rec['Median'] = '-'
        del rec['folds']

        rec['Samples'] = rec['S'] + rec['I'] + rec['R']

    records = list(records.values())

    records = apply_modifier(records, record_modifier)

    return records


def apply_modifier(records, record_modifier):
    results = []
    if record_modifier:
        for rec in records:
            results.append(
                record_modifier(rec)
            )
    else:
        results = records

    return results


def record_modifier(record):
    variant_name = record['Variant name']
    reference = record['Reference']

    if variant_name.endswith('full genome'):
        reference = '{}*'.format(reference)
        variant_name = variant_name.split()[0]

    record['Variant name'] = variant_name
    record['Reference'] = reference
    return record


def gen_plasma_aggre_table(
        conn, row_filters, subrow_filters,
        sql_template, record_modifier=None):

    cursor = conn.cursor()

    records = []
    for row_name, attr_r in row_filters.items():
        for subrow_name, attr_subr in subrow_filters.items():
            rxtype = attr_subr['rxtype']

            r_filter = attr_r.get('filter', [])
            filter = '\n    '.join(r_filter)

            if subrow_name.lower().startswith('cp'):
                filter += '\n   '
                filter += '\n   '.join(attr_subr.get('cp_filters', []))

            sql = sql_template.format(
                rxtype=rxtype,
                filters=filter
            )
            # print(sql)

            cursor.execute(sql)

            groups = defaultdict(list)
            for row in cursor.fetchall():
                variant_name = row_name
                control = row['control']
                cp_name = row['rx_name']
                reference = row['ref_name']

                cp_name = rename_cp(cp_name, reference)

                group_key = '{}{}{}{}'.format(
                    variant_name,
                    control,
                    cp_name,
                    reference
                )
                groups[group_key].append(row)

            for _, r_list in groups.items():
                cp_name = r_list[0]['rx_name']
                reference = r_list[0]['ref_name']
                cp_name = rename_cp(cp_name, reference)

                all_fold = [
                    [r['fold']] * r['sample_count']
                    for r in r_list if r['fold']]
                all_fold = [r for j in all_fold for r in j]
                s_fold = [r for r in all_fold if is_susc(r)]
                i_fold = [
                    r for r in all_fold if is_partial_resistant(r)]
                r_fold = [r for r in all_fold if is_resistant(r)]

                num_s_fold = len(s_fold)
                num_i_fold = len(i_fold)
                num_r_fold = len(r_fold)
                num_results = num_s_fold + num_i_fold + num_r_fold
                median_fold = median(all_fold)

                rec = {
                    'Variant name': variant_name,
                    'Plasma': PLASMA_POST_RENAME.get(cp_name, cp_name),
                    'Samples': num_results,
                    'Reference': reference,
                    'Median': median_fold,
                    'S': num_s_fold,
                    'I': num_i_fold,
                    'R': num_r_fold,
                }

                records.append(rec)

    records = apply_modifier(records, record_modifier)

    return records


def convert_to_json(json_save_path, records):
    json_results = defaultdict(list)
    for r in records:
        variant = r['Variant name']
        json_results[variant].append({
            'variant': variant,
            'rx': r['Plasma'],
            'samples': r['Samples'],
            's_fold': r['S'],
            'i_fold': r['I'],
            'r_fold': r['R'],
            'reference': r['Reference'],
            'median': r['Median'],
        })

    results = []
    for variant, assays in json_results.items():
        results.append({
            'variant': variant,
            'assays': sorted(assays, key=itemgetter('rx')),
        })

    variant = sorted(results, key=itemgetter('variant'))
    dump_json(json_save_path, results)


def rename_cp(cp_name, reference):
    cp_name = PLASMA_RENAME.get(cp_name, cp_name)
    rename_executors = RENAME_CP_EXECUTOR.get(reference, [])
    for tester, new_name in rename_executors:
        if tester(cp_name):
            cp_name = new_name

    return cp_name
