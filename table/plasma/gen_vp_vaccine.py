from preset import dump_csv
from preset import dump_json
from preset import DATA_FILE_PATH
from preset import round_number
from preset import group_records_by
from sql import row2dict
from operator import itemgetter
from statistics import median
from resistancy import is_susc
from resistancy import is_partial_resistant
from resistancy import is_resistant

SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.fold,
    rx.vaccine_name,
    rx.dosage,
    rx.infected_var_name,
    SUM(s.cumulative_count) num_fold,
    vac.priority
FROM
    susc_results_view s,
    rx_vacc_plasma rx,
    vaccines vac
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.fold IS NOT NULL
    AND
    rx.vaccine_name = vac.vaccine_name
GROUP BY
    s.ref_name,
    s.rx_name,
    s.control_iso_name,
    s.iso_name
"""


def gen_vp_vaccine(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)
    records = row2dict(cursor.fetchall())

    [
        r.update({
            'infection': (
                'infected' if r['infected_var_name']
                else '')
        })
        for r in records
    ]

    num_records = sum([r['num_fold'] for r in records])

    vaccine_group = group_records_by(records, 'vaccine_name')
    vaccine_dosage_group = {}
    for vaccine, rx_list in vaccine_group.items():
        for dosage, dosage_list in group_records_by(rx_list, 'dosage').items():
            vaccine_dosage_group[(vaccine, dosage)] = dosage_list

    vaccine_dosage_infected = {}
    for (vaccine, dosage), dosage_list in vaccine_dosage_group.items():
        for infected, infected_list in group_records_by(
                dosage_list, 'infection').items():
            vaccine_dosage_infected[
                (vaccine, dosage, infected)] = infected_list

    vaccine_results = []
    for (vaccine, dosage, infected), rx_list in vaccine_dosage_infected.items():

        all_fold = [[r['fold']] * r['num_fold'] for r in rx_list]
        all_fold = [r for j in all_fold for r in j]
        num_fold = len(all_fold)

        s_fold = [r for r in all_fold if is_susc(r)]
        i_fold = [r for r in all_fold if is_partial_resistant(r)]
        r_fold = [r for r in all_fold if is_resistant(r)]

        num_s_fold = len(s_fold)
        num_i_fold = len(i_fold)
        num_r_fold = len(r_fold)
        median_fold = round_number(median(all_fold))

        vaccine_results.append({
            'vaccine': vaccine,
            'dosage': dosage,
            'infection': infected,
            'num_ref_name': len(set(
                r['ref_name'] for r in rx_list
            )),
            'num_fold': num_fold,
            'S': num_s_fold,
            'I': num_i_fold,
            'R': num_r_fold,
            'median_fold': float(median_fold),
            'percent': float(
                round_number(num_fold / num_records * 100)),
            'priority': rx_list[0]['priority']
        })
    vaccine_results.sort(
        key=itemgetter('priority', 'dosage', 'infection'))
    vaccine_results = [
        i for i in vaccine_results
        if i['num_ref_name'] != 1 and (
            i['vaccine'] not in ['Inactivated', 'mRNA']
        )
    ]
    dump_csv(DATA_FILE_PATH / 'table_vp_vaccine.csv', vaccine_results)
    dump_json(DATA_FILE_PATH / 'table_vp_vaccine.json', vaccine_results)
