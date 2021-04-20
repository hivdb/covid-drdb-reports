from collections import defaultdict
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import MAB_RENAME
from mab.preset import RX_MAB
from variant.preset import CONTROL_VARIANTS_SQL
from variant.preset import filter_by_variant

SQL = """
SELECT
    rx.ab_name,
    s.cumulative_count as count,
    s.variant_name,
    rx.availability as avail,
    rx.pdb_id as pdb,
    rx.target as target
FROM
    susc_results as s,
    ({rx_type}) as rx
ON
    s.ref_name = rx.ref_name AND
    s.rx_name = rx.rx_name
WHERE
    s.control_variant_name in {control_variants}
    AND s.fold IS NOT NULL
""".format(rx_type=RX_MAB, control_variants=CONTROL_VARIANTS_SQL)


def gen_table_all_mab(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)

    mab_group = defaultdict(list)
    records = cursor.fetchall()
    records = filter_by_variant(records)
    for rec in records:
        ab_name = rec['ab_name']
        ab_name = MAB_RENAME.get(ab_name, ab_name)

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

    save_path = DATA_FILE_PATH / 'summary_mab.csv'
    dump_csv(save_path, record_list)
