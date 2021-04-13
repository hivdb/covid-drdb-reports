from collections import defaultdict
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import MAB_RENAME
from mab.preset import RX_MAB

SQL = """
SELECT
    rx.ab_name,
    s.cumulative_count as count,
    rx.availability as avail,
    rx.pdb_id as pdb,
    rx.target as target
FROM
    susc_results as s,
    ({rx_type}) as rx
ON
    s.ref_name = rx.ref_name AND
    s.rx_name = rx.rx_name
""".format(rx_type=RX_MAB)


def gen_table_all_mab(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)

    mab_group = defaultdict(list)
    for rec in cursor.fetchall():
        ab_name = rec['ab_name']
        ab_name = MAB_RENAME.get(ab_name, ab_name)

        mab_group[ab_name].append(rec)

    record_list = []
    for ab_name, rlist in mab_group.items():
        count = sum([r['count'] for r in rlist] + [0])
        target = rlist[0]['target']
        pdb = rlist[0]['pdb']
        avail = rlist[0]['avail']

        record_list.append({
            'mab': ab_name,
            'results': count,
            'avail': avail or '',
            'target': target or '',
            'pdb': pdb or '',
        })

    record_list.sort(key=itemgetter(
        'avail',
        'target',
        'results',
        'pdb',
        ), reverse=True)

    save_path = DATA_FILE_PATH / 'table_mab_figure.csv'
    dump_csv(save_path, record_list)
