from collections import defaultdict
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import get_fold_stat


SQL_TMPL = """
SELECT
    s.iso_name,
    s.fold_cmp,
    s.fold,
    s.cumulative_count num_fold,
    rx.ab_name,
    rx.availability avail,
    rx.pdb_id pdb,
    rx.target target,
    iso.*
FROM
    susc_results s,
    rx_mab_view rx,
    {iso_type} iso
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.iso_name = iso.iso_name
    ;
"""


def gen_table_variant_mab(conn):
    iso_type = 'isolate_mutations_single_s_mut_view'
    save_path = DATA_FILE_PATH / 'variant' / 'summary_single_mab.csv'
    by_single(conn, iso_type, save_path)

    iso_type = 'isolate_mutations_combo_s_mut_view'
    save_path = DATA_FILE_PATH / 'variant' / 'summary_combo_mab.csv'
    by_combo(conn, iso_type, save_path)


def by_single(conn, iso_type, save_path):
    sql = SQL_TMPL.format(iso_type=iso_type)
    cursor = conn.cursor()
    cursor.execute(sql)
    db_records = cursor.fetchall()

    mut_group = defaultdict(list)
    for rec in db_records:
        mut_name = rec['single_mut_name']
        mut_group[mut_name].append(rec)

    record_list = []
    for mut_name, rlist in mut_group.items():
        rx_group = defaultdict(list)
        for rec in rlist:
            ab_name = rec['ab_name']
            rx_group[ab_name].append(rec)

        for ab_name, rx_list in rx_group.items():
            avail = rx_list[0]['avail']
            target = rx_list[0]['target']

            num_s, num_i, num_r, median_fold, num_fold = get_fold_stat(rx_list)

            record_list.append({
                'pattern': mut_name,
                'ref': rx_list[0]['ref'],
                'pos': rx_list[0]['position'],
                'aa': rx_list[0]['amino_acid'],
                'domain': rx_list[0]['domain'],
                'mab': ab_name,
                'avail': avail,
                'target': target,
                'median_fold': median_fold,
                'num_fold': num_fold,
                'S': num_s,
                'I': num_i,
                'R': num_r,
            })

    record_list.sort(key=itemgetter(
        'pos',
        'aa',
        'mab',
        'avail',
        ))
    dump_csv(save_path, record_list)


def by_combo(conn, iso_type, save_path):
    sql = SQL_TMPL.format(iso_type=iso_type)
    cursor = conn.cursor()
    cursor.execute(sql)
    db_records = cursor.fetchall()

    mut_group = defaultdict(list)
    for rec in db_records:
        pattern = rec['pattern']
        mut_group[pattern].append(rec)

    record_list = []
    for pattern, rlist in mut_group.items():
        rx_group = defaultdict(list)
        for rec in rlist:
            ab_name = rec['ab_name']
            rx_group[ab_name].append(rec)

        for ab_name, rx_list in rx_group.items():
            avail = rx_list[0]['avail']
            target = rx_list[0]['target']

            num_s, num_i, num_r, median_fold, num_fold = get_fold_stat(rx_list)

            record_list.append({
                'pattern': pattern,
                'var_name': rx_list[0]['var_name'] or '',
                'mab': ab_name,
                'avail': avail,
                'target': target,
                'median_fold': median_fold,
                'num_fold': num_fold,
                'S': num_s,
                'I': num_i,
                'R': num_r,
            })

    record_list.sort(key=itemgetter(
        'var_name',
        'pattern',
        'mab',
        'avail',
        ))

    dump_csv(save_path, record_list)
