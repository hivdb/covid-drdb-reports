from preset import DATA_FILE_PATH
from preset import dump_csv
from collections import defaultdict
from sql import row2dict

SQL_TMPL = """
SELECT
    s.ref_name ref_name,
    s.control_iso_name control,
    s.iso_name iso_name,
    s.cumulative_count num_fold,
    s.fold_cmp fold_cmp,
    s.fold fold,
    iso.*
FROM
    susc_results_aggr_50_wt_view s,
    {rx_type} rx,
    {iso_type} iso
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.iso_name = iso.iso_name
;
"""


def gen_table_variant_aggre(conn):
    rx_type = 'rx_vacc_plasma'
    save_path = DATA_FILE_PATH / 'variant' / 'summary_vp_aggre.csv'
    _gen_table_variant_aggre(conn, rx_type, save_path)

    rx_type = 'rx_conv_plasma'
    save_path = DATA_FILE_PATH / 'variant' / 'summary_cp_aggre.csv'
    _gen_table_variant_aggre(conn, rx_type, save_path)


def _gen_table_variant_aggre(conn, rx_type, save_path):
    results = []

    cursor = conn.cursor()
    iso_type = 'isolate_mutations_single_s_mut_view'
    sql = SQL_TMPL.format(
        rx_type=rx_type,
        iso_type=iso_type
        )

    # print(sql)
    cursor.execute(sql)
    records = row2dict(cursor.fetchall())
    results += get_fold_by_single(records)

    iso_type = 'isolate_mutations_combo_s_mut_view'
    sql = SQL_TMPL.format(
        rx_type=rx_type,
        iso_type=iso_type
        )

    # print(sql)
    cursor.execute(sql)
    records = row2dict(cursor.fetchall())
    results += get_fold_by_combo(records)

    dump_csv(save_path, results)


def get_fold_by_single(rx_list):
    mut_group = defaultdict(list)
    for rec in rx_list:
        mut_name = rec['single_mut_name']
        mut_group[mut_name].append(rec)

    results = []
    for mut_name, rx_list in mut_group.items():
        for r in rx_list:
            results.append({
                'pattern': mut_name,
                'domain': r.get('domain'),
                'var_name': r.get('var_name'),
                'ref_name': r['ref_name'],
                'fold_cmp': r['fold_cmp'],
                'median': r['fold'],
                'num_fold': r['num_fold']
            })

    return results


def get_fold_by_combo(rx_list):

    mut_group = defaultdict(list)
    for rec in rx_list:
        mut_name = rec['pattern']
        mut_group[mut_name].append(rec)

    results = []
    for mut_name, r_list in mut_group.items():
        for r in r_list:
            results.append({
                'pattern': mut_name,
                'domain': r.get('domain'),
                'var_name': r.get('var_name'),
                'ref_name': r['ref_name'],
                'fold_cmp': r['fold_cmp'],
                'median': r['fold'],
                'num_fold': r['num_fold']
            })

    return results
