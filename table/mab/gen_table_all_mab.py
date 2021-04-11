from collections import defaultdict
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import MAB_RENAME
from mab.preset import ANTIBODY_TARGET_SQL

SQL = """
SELECT
    s.rx_name,
    s.cumulative_count as count,
    a.availability as avail,
    t.pdb_id as pdb,
    t.target as target
FROM
    susc_results as s
INNER JOIN rx_antibodies as r ON
    s.ref_name = r.ref_name
    AND s.rx_name = r.rx_name
LEFT JOIN antibodies as a ON a.ab_name = r.ab_name
LEFT JOIN """ + ANTIBODY_TARGET_SQL + """ as t ON a.ab_name = t.ab_name
WHERE
    r.ab_name in (
        SELECT ab_name FROM 'antibodies'
    )
"""


def gen_table_all_mab(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)

    mab_group = defaultdict(list)
    for rec in cursor.fetchall():
        rx_name = rec['rx_name']
        rx_name = MAB_RENAME.get(rx_name, rx_name)

        mab_group[rx_name].append(rec)

    record_list = []
    for rx_name, rlist in mab_group.items():
        count = sum([r['count'] for r in rlist] + [0])
        target = rlist[0]['target']
        pdb = rlist[0]['pdb']
        avail = rlist[0]['avail']

        record_list.append({
            'mab': rx_name,
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
