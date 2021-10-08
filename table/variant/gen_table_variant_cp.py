from collections import defaultdict
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import get_fold_stat, group_var_name

SQL_TMPL = """
SELECT
    s.iso_name,
    s.ref_name,
    s.fold_cmp,
    s.fold,
    s.cumulative_count as num_fold,
    iso.*
FROM
    susc_results s,
    rx_conv_plasma r,
    {iso_type} iso
WHERE
    s.ref_name = r.ref_name
    AND
    s.rx_name = r.rx_name
    AND
    s.iso_name = iso.iso_name
;
"""


def gen_table_variant_cp(conn):
    iso_type = 'isolate_mutations_single_s_mut_view'
    save_path = DATA_FILE_PATH / 'variant' / 'summary_single_cp.csv'
    by_single(conn, iso_type, save_path)

    iso_type = 'isolate_mutations_combo_s_mut_view'
    save_path = DATA_FILE_PATH / 'variant' / 'summary_combo_cp.csv'
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
    for mut_name, rx_list in mut_group.items():
        num_s, num_i, num_r, median_fold, num_fold = get_fold_stat(rx_list)

        record_list.append({
            'pattern': mut_name,
            'ref': rx_list[0]['ref'],
            'pos': rx_list[0]['position'],
            'aa': rx_list[0]['amino_acid'],
            'domain': rx_list[0]['domain'],
            'median_fold': median_fold,
            'num_ref_name': len(set([
                r['ref_name']
                for r in rx_list
            ])),
            'num_fold': num_fold,
            'S': num_s,
            'I': num_i,
            'R': num_r,
        })

    record_list.sort(key=itemgetter(
        'pos',
        'aa',
        ))

    record_list.append({
        'pattern': 'summary',
        'num_ref_name': len(set([
            r['ref_name']
            for r in db_records
        ])),
        'num_fold': sum([r['num_fold'] for r in db_records] + [0]),
    })

    dump_csv(save_path, record_list)


def by_combo(conn, iso_type, save_path):
    sql = SQL_TMPL.format(iso_type=iso_type)
    cursor = conn.cursor()
    cursor.execute(sql)
    db_records = cursor.fetchall()

    mut_group = defaultdict(list)
    for rec in db_records:
        var_name = rec['var_name']
        var_name = group_var_name(var_name)
        mut_group[var_name].append(rec)

    record_list = []
    for var_name, rx_list in mut_group.items():
        num_s, num_i, num_r, median_fold, num_fold = get_fold_stat(rx_list)

        record_list.append({
            'pattern': var_name,
            'var_name': var_name,
            'median_fold': median_fold,
            'num_ref_name': len(set([
                r['ref_name']
                for r in rx_list
            ])),
            'num_fold': num_fold,
            'S': num_s,
            'I': num_i,
            'R': num_r,
        })

    record_list.sort(key=itemgetter(
        'var_name',
        'pattern',
        ))

    record_list.append({
        'pattern': 'summary',
        'num_ref_name': len(set([
            r['ref_name']
            for r in db_records
        ])),
        'num_fold': sum([r['num_fold'] for r in db_records] + [0]),
    })

    dump_csv(save_path, record_list)


