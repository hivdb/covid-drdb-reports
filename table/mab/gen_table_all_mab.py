from collections import defaultdict
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import MAB_RENAME
from mab.preset import RX_MAB
from variant.preset import filter_by_variant

SQL = """
SELECT
    rx.ab_name,
    rx.synonyms,
    s.cumulative_count AS count,
    s.iso_name,
    rx.availability AS avail,
    rx.pdb_id AS pdb,
    rx.target AS target,
    rx.class AS class,
    rx.epitope AS epitope,
    iso.var_name
FROM
    susc_results AS s,
    ({rx_type}) AS rx,
    isolates AS iso
ON
    s.ref_name = rx.ref_name AND
    s.rx_name = rx.rx_name AND
    s.iso_name = iso.iso_name
WHERE
    s.potency_type LIKE 'IC%'
;
""".format(rx_type=RX_MAB)


def gen_table_all_mab(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)

    mab_group = defaultdict(list)
    records = cursor.fetchall()
    records = filter_by_variant(records)

    for rec in records:
        ab_name = rec['ab_name']
        ab_name = MAB_RENAME.get(ab_name, ab_name)
        mab_group[ab_name].append(rec)

    record_list = []
    for ab_name, rlist in mab_group.items():
        count = sum([r['count'] for r in rlist] + [0])
        target = rlist[0]['target']
        pdb = rlist[0]['pdb']
        avail = rlist[0]['avail']
        synonyms = rlist[0]['synonyms']
        epitope = rlist[0]['epitope']
        ab_class = rlist[0]['class']

        record_list.append({
            'mab': (
                ab_name if not ab_name[0].isdigit() else "'{}".format(ab_name)
            ),
            'synonyms': synonyms,
            'results': count,
            'avail': avail or '',
            'target': target or '',
            'class': ab_class or '',
            'pdb': pdb or '',
            'epitope': epitope or '',
        })

    record_list.sort(key=itemgetter(
        'avail',
        'target',
        'results',
        'pdb',
        ), reverse=True)

    save_path = DATA_FILE_PATH / 'summary_mab.csv'
    dump_csv(save_path, record_list)
