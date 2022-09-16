from sql import get_sql_stmt
from pathlib import Path
from sql import row2dict
from preset import group_records_by
from preset import dump_csv, load_csv
from scipy.stats.mstats import gmean
from preset import DATA_FILE_PATH


def analyze_cp_titer(conn):

    table = load_csv(
        DATA_FILE_PATH / 'figure' / 'figure_plasma_variant_titer.csv')

    table = [
        i
        for i in table
        if not i['dosage']
    ]

    report = []

    for inf, inf_rec_list in group_records_by(
            table, 'infections').items():
        for omic, omic_rec_list in group_records_by(
                inf_rec_list, 'var_name').items():
            for month, month_rec in group_records_by(
                    omic_rec_list, 'month').items():
                for ref, ref_list in group_records_by(
                        month_rec, 'ref_name').items():
                    report.append({
                        'infection': inf,
                        'var_name': omic,
                        'month': month,
                        'ref': ref,
                        'num_titer': sum([
                            i['num_exp'] for i in ref_list]),
                        'titer': gmean([
                            i['titer'] for i in ref_list],
                            weights=[
                                i['num_exp'] for i in ref_list])
                    })

    save_file_path = DATA_FILE_PATH / 'omicron_plasma' / 'cp_titer.csv'
    dump_csv(save_file_path, report)

    _stat_analysis(save_file_path)


def _stat_analysis(file_path):

    table = load_csv(file_path)

    report = []
    for inf, inf_rec_list in group_records_by(
            table, 'infection').items():
        for omic, omic_rec_list in group_records_by(
                inf_rec_list, 'var_name').items():
            for month, month_rec in group_records_by(
                    omic_rec_list, 'month').items():
                m = gmean(
                    [i['titer'] for i in month_rec],
                    weights=[
                        i['num_titer'] for i in month_rec])

                # st = stdev([
                #     i['titer']
                #     for i in month_rec
                # ]) if len(month_rec) > 1 else ''

                for rec in month_rec:
                    rec['mean'] = m
                    # rec['z-score'] = (rec['titer'] - m) / st if st else ''
                    # if (
                    #     rec['z-score'] and (
                    #         rec['z-score'] < -2 or rec['z-score'] > 2)):
                    #     rec['outlier'] = True
                    report.append(rec)

    dump_csv(DATA_FILE_PATH / 'omicron_plasma' / 'cp_titer_stat.csv', report)
