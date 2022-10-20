from collections import defaultdict
from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json
from operator import itemgetter

SQL = """
SELECT
    rx.ab_name,
    rx.synonyms,
    s.ref_name,
    s.cumulative_count num_fold,
    rx.availability avail,
    rx.pdb_id pdb,
    rx.target target,
    rx.class class,
    rx.epitope epitope,
    rx.institute institute,
    rx.origin origin,
    iso.var_name
FROM
    susc_results_view s,
    rx_mab_view rx,
    isolates iso
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.iso_name = iso.iso_name
;
"""


def mab_summary(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)

    mab_group = defaultdict(list)
    records = cursor.fetchall()

    for rec in records:
        ab_name = rec['ab_name']
        mab_group[ab_name].append(rec)

    record_list = []
    for ab_name, rlist in mab_group.items():
        record_list.append({
            'mab': (
                ab_name if not ab_name[0].isdigit() else "'{}".format(ab_name)
            ),
            'synonyms': rlist[0]['synonyms'],
            'num_ref_name': len(set(
                r['ref_name'] for r in rlist
            )),
            'num_fold': sum([r['num_fold'] for r in rlist] + [0]),
            'avail': rlist[0]['avail'] or '',
            'target': rlist[0]['target'] or '',
            'class': rlist[0]['class'] or '',
            'pdb': rlist[0]['pdb'] or '',
            'institute': rlist[0]['institute'] or '',
            'origin': rlist[0]['origin'] or '',
            # 'epitope': rlist[0]['epitope'] or '',
        })

    record_list.sort(key=itemgetter(
        'avail',
        'num_fold',
        'mab',
        ), reverse=True)

    save_path = DATA_FILE_PATH / 'mab' / 'summary_mab.csv'
    dump_csv(save_path, record_list)

    record_list = [
        i
        for i in record_list
        if i['synonyms'] and i['avail']
    ]

    dump_json(DATA_FILE_PATH / 'table_mab_detail.json', record_list)
