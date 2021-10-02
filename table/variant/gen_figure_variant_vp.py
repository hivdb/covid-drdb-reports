from collections import defaultdict
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter

SQL_TMPL = """
SELECT
    susc.fold_cmp,
    susc.fold,
    susc.cumulative_count as count,
    susc.vaccine_name,
    susc.vaccine_type,
    iso.*
FROM
    susc_results_vp_50_wt_view susc,
    {iso_type} iso
WHERE
    susc.iso_name = iso.iso_name
"""


def gen_figure_variant_vp(conn):
    iso_type = 'isolate_mutations_single_s_mut_view'
    save_path = DATA_FILE_PATH / 'figure' / 'variant_single_vp.csv'
    by_single(conn, iso_type, save_path)

    iso_type = 'isolate_mutations_combo_s_mut_view'
    save_path = DATA_FILE_PATH / 'figure' / 'variant_combo_vp.csv'
    by_combo(conn, iso_type, save_path)


def by_single(conn, iso_type, save_path):
    sql = SQL_TMPL.format(iso_type=iso_type)

    cursor = conn.cursor()

    cursor.execute(sql)

    results = []
    for rec in cursor.fetchall():
        results.append({
            'pattern': rec['single_mut_name'],
            'vaccine': rec['vaccine_name'],
            'vaccine_type': rec['vaccine_type'],
            'ref': rec['ref'],
            'pos': rec['position'],
            'aa': rec['amino_acid'],
            'domain': rec['domain'],
            'fold': rec['fold'],
        })

    results.sort(key=itemgetter('pos', 'aa', 'vaccine'))
    dump_csv(save_path, results)


def by_combo(conn, iso_type, save_path):
    sql = SQL_TMPL.format(iso_type=iso_type)

    cursor = conn.cursor()

    cursor.execute(sql)

    results = []
    for rec in cursor.fetchall():
        results.append({
            'pattern': rec['pattern'],
            'var_name': rec['var_name'] or '',
            'vaccine': rec['vaccine_name'],
            'vaccine_type': rec['vaccine_type'],
            'fold': rec['fold'],
        })

    results.sort(key=itemgetter('var_name', 'pattern', 'vaccine'))
    dump_csv(save_path, results)
