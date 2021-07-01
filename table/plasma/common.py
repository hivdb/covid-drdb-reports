import statistics
from operator import itemgetter
from collections import defaultdict
from resistancy import round_fold
from preset import dump_json
from resistancy import RESISTANCE_FILTER
from resistancy import is_susc
from resistancy import is_partial_resistant
from resistancy import is_resistant
from resistancy import get_susceptibility
from statistics import median
from variant.preset import CONTROL_VARIANTS_SQL


def gen_plasma_indiv_table(
        conn, row_filters, subrow_filters,
        sql_template,
        plasma_type=None, record_modifier=None):

    cursor = conn.cursor()

    records = defaultdict(dict)
    for iso_name, attr_r in row_filters.items():
        for plasma_name, attr_subr in subrow_filters.items():
            for resist_name, resist_filter in RESISTANCE_FILTER.items():
                rxtype = attr_subr['rxtype']

                r_filter = attr_r.get('filter', [])
                filter = '\n    '.join(r_filter + resist_filter)

                if plasma_name.lower().startswith('cp'):
                    filter += '\n   '
                    filter += '\n   '.join(attr_subr.get('cp_filters', []))

                sql = sql_template.format(
                    rxtype=rxtype,
                    filters=filter,
                    control_variants=CONTROL_VARIANTS_SQL,
                )
                # print(sql)

                cursor.execute(sql)
                for row in cursor.fetchall():
                    cp_name = row['rx_name']
                    reference = row['ref_name']
                    num_results = row['sample_count']
                    fold = row['fold']

                    key = '{}{}{}'.format(iso_name, cp_name, reference)

                    if plasma_type == 'VP':
                        dosage = row['dosage']
                        key = '{}{}{}{}'.format(
                            iso_name, cp_name, reference, dosage)

                    rec = records[key]
                    rec['Variant name'] = iso_name
                    rec['Plasma'] = cp_name
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
                    rec['Aggregate'] = False

                    if plasma_type == 'VP':
                        rec['dosage'] = row['dosage']

    for rec in records.values():
        folds = rec.get('folds', [])
        if folds:
            rec['Median'] = str(round_fold(statistics.median(folds)))
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
    iso_name = record['Variant name']
    reference = record['Reference']

    if iso_name.endswith('full genome'):
        reference = '{}*'.format(reference)
        iso_name = iso_name.split()[0]

    record['Variant name'] = iso_name
    record['Reference'] = reference
    return record


def gen_plasma_aggre_table(
        conn, row_filters, subrow_filters,
        sql_template,
        plasma_type=None, record_modifier=None):

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
                filters=filter,
                control_variants=CONTROL_VARIANTS_SQL,
            )
            # print(sql)

            cursor.execute(sql)

            groups = defaultdict(list)
            for row in cursor.fetchall():
                iso_name = row_name
                control = row['control']
                cp_name = row['rx_name']
                reference = row['ref_name']

                group_key = '{}{}{}{}'.format(
                    iso_name,
                    control,
                    cp_name,
                    reference
                )
                if plasma_type == 'VP':
                    dosage = row['dosage']
                    group_key = '{}{}{}{}{}'.format(
                        iso_name,
                        control,
                        cp_name,
                        reference,
                        dosage
                    )
                groups[group_key].append(row)

            for _, r_list in groups.items():
                cp_name = r_list[0]['rx_name']
                reference = r_list[0]['ref_name']

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
                    'Variant name': iso_name,
                    'Plasma': cp_name,
                    'Samples': num_results,
                    'Reference': reference,
                    'Median': median_fold,
                    'S': num_s_fold,
                    'I': num_i_fold,
                    'R': num_r_fold,
                    'Aggregate': True,
                }

                if plasma_type == 'VP':
                    rec['dosage'] = r_list[0]['dosage']

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


def get_sample_number_pair(indiv_list, aggre_list):

    num_s = sum([int(r['S']) for r in indiv_list])
    num_i = sum([int(r['I']) for r in indiv_list])
    num_r = sum([int(r['R']) for r in indiv_list])

    aggre_s = [
        r for r in aggre_list
        if get_susceptibility(r['Fold']) == 'susceptible']
    aggre_i = [
        r for r in aggre_list
        if get_susceptibility(r['Fold']) == 'partial-resistance']
    aggre_r = [
        r for r in aggre_list
        if get_susceptibility(r['Fold']) == 'resistant']

    num_s_aggre = sum([int(r['Samples']) for r in aggre_s])
    num_i_aggre = sum([int(r['Samples']) for r in aggre_i])
    num_r_aggre = sum([int(r['Samples']) for r in aggre_r])

    return (num_s + num_s_aggre), (num_i + num_i_aggre), (num_r + num_r_aggre)
