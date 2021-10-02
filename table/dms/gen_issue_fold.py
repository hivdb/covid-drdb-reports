from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import load_csv


def gen_issue_fold(
        input_file_path=DATA_FILE_PATH / 'dms' / 'fold_dms.csv',
        output_file_path=DATA_FILE_PATH / 'dms' / 'fold_dms_issue.csv'):

    results = []

    for rec in load_csv(input_file_path):
        fold = float(rec['fold'])
        score = float(rec['score'])

        if fold <= 10 and score <= 0.1:
            continue
        if fold > 10 and score > 0.1:
            continue

        results.append(rec)

    dump_csv(output_file_path, results)
