from collections import defaultdict
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter


SQL_TMPL = """
SELECT
    susc.ab_name,
    susc.fold_cmp,
    susc.fold,
    susc.cumulative_count as count,
    susc.availability as avail,
    susc.pdb_id as pdb,
    susc.target as target,
    iso.*
FROM
    susc_results_mab_50_wt_view as susc,
    {iso_type} iso
WHERE
    susc.iso_name = iso.iso_name
"""


def gen_figure_variant_mab(conn):

    iso_type = 'isolate_mutations_single_s_mut_view'
    save_path = DATA_FILE_PATH / 'figure' / 'variant_single_mab.csv'
    by_single(conn, iso_type, save_path)

    iso_type = 'isolate_mutations_combo_s_mut_view'
    save_path = DATA_FILE_PATH / 'figure' / 'variant_combo_mab.csv'
    by_combo(conn, iso_type, save_path)


def by_single(conn, iso_type, save_path):
    sql = SQL_TMPL.format(iso_type=iso_type)

    cursor = conn.cursor()

    cursor.execute(sql)

    results = []
    for rec in cursor.fetchall():
        results.append({
            'pattern': rec['single_mut_name'],
            'ref': rec['ref'],
            'pos': rec['position'],
            'aa': rec['amino_acid'],
            'domain': rec['domain'],
            'mab': rec['ab_name'],
            'avail': rec['avail'],
            'target': rec['target'],
            'fold': rec['fold'],
        })

    results.sort(key=itemgetter(
        'pos',
        'aa',
        'mab',
        'avail'))
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
            'mab': rec['ab_name'],
            'avail': rec['avail'],
            'target': rec['target'],
            'fold': rec['fold'],
        })

    results.sort(key=itemgetter(
        'var_name',
        'pattern',
        'mab',
        'avail',))
    dump_csv(save_path, results)
