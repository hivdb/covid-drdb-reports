from collections import defaultdict
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import MAB_RENAME

SQL = """
SELECT
    s.rx_name,
    s.cumulative_count as count
FROM
    susc_results as s,
    rx_antibodies as r
ON
    s.ref_name = r.ref_name
    AND s.rx_name = r.rx_name
WHERE
    r.ab_name in (
        SELECT ab_name FROM 'antibodies'
    )
"""


# 3951

def gen_table_all_mab(conn):
    cursor = conn.cursor()

    cursor.execute(SQL)

    mab_group = defaultdict(list)
    for i in cursor.fetchall():
        rx_name = i['rx_name']
        rx_name = MAB_RENAME.get(rx_name, rx_name)
        num = i['count']
        mab_group[rx_name].append({
            'rx_name': rx_name,
            'num': num
        })

    record_list = []
    for rx_name, rlist in mab_group.items():
        count = 1
        for rec in rlist:
            count += rec['num']

        record_list.append({
            'mab': rx_name,
            'count': count
        })

    record_list.sort(key=itemgetter('count'), reverse=True)

    save_path = DATA_FILE_PATH / 'table_mab_figure.csv'
    dump_csv(save_path, record_list)
