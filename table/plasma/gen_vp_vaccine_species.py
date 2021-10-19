from preset import dump_csv
from preset import DATA_FILE_PATH
from collections import defaultdict
from operator import itemgetter

SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    rx.vaccine_name,
    sub.subject_species species,
    SUM(s.cumulative_count) num_fold
FROM
    susc_results_view s,
    rx_vacc_plasma rx,
    subjects sub
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    rx.ref_name = sub.ref_name
    AND
    rx.subject_name = sub.subject_name
    AND
    s.fold IS NOT NULL
GROUP BY
    s.ref_name,
    s.rx_name,
    s.control_iso_name,
    s.iso_name
"""


def gen_vp_vaccine_species(conn):
    cursor = conn.cursor()
    cursor.execute(SQL)
    records = cursor.fetchall()

    vaccine_species_group = defaultdict(list)
    for rec in records:
        vaccine = rec['vaccine_name']
        species = rec['species']
        vaccine_species_group[(vaccine, species)].append(rec)

    vaccine_species_results = []
    for (vaccine, species), rx_list in vaccine_species_group.items():
        vaccine_species_results.append({
            'vaccine': vaccine,
            'num_ref_name': len(set(
                r['ref_name'] for r in rx_list
            )),
            'num_fold': sum([r['num_fold'] for r in rx_list]),
            'species': species,
        })
    vaccine_species_results.sort(key=itemgetter('vaccine', 'species'))
    save_path = DATA_FILE_PATH / 'vp' / 'summary_vp_vaccine_species.csv'
    dump_csv(save_path, vaccine_species_results)
