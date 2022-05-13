from collections import defaultdict
from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json
from operator import itemgetter
from .preset import get_fold_stat, group_var_name


SQL_TMPL = """
SELECT
    s.iso_name,
    s.ref_name,
    s.fold_cmp,
    s.fold,
    s.cumulative_count num_fold,
    r.vaccine_name,
    vac.vaccine_type,
    iso.*
FROM
    susc_results_view s,
    rx_vacc_plasma r,
    vaccines vac,
    {iso_type} iso
WHERE
    s.ref_name = r.ref_name
    AND
    s.rx_name = r.rx_name
    AND
    r.vaccine_name = vac.vaccine_name
    AND
    s.iso_name = iso.iso_name
;
"""


def gen_table_variant_vp(conn):
    iso_type = 'isolate_mutations_single_s_mut_view'
    save_path = DATA_FILE_PATH / 'variant' / 'summary_single_vp.csv'
    single_mut_records = by_single(conn, iso_type, save_path)

    iso_type = 'isolate_mutations_combo_s_mut_view'
    save_path = DATA_FILE_PATH / 'table_vp_variant.csv'
    by_combo(conn, iso_type, save_path, single_mut_records)


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
        # vacc_group = defaultdict(list)
        # for rec in rx_list:
        #     vacc_name = rec['vaccine_name']
        #     vacc_group[vacc_name].append(rec)

        # for vacc_name, rx_list in vacc_group.items():
        num_s, num_i, num_r, median_fold, num_fold = get_fold_stat(rx_list)

        record_list.append({
            'pattern': mut_name,
            # 'vaccine': rx_list[0]['vaccine_name'],
            # 'vaccine_type': rx_list[0]['vaccine_type'],
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
        # 'vaccine',
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

    return db_records


def by_combo(conn, iso_type, save_path, single_mut_records):
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
        # vacc_group = defaultdict(list)
        # for rec in rx_list:
        #     vacc_name = rec['vaccine_name']
        #     vacc_group[vacc_name].append(rec)

        # for vacc_name, rx_list in vacc_group.items():
        num_s, num_i, num_r, median_fold, num_fold = get_fold_stat(rx_list)

        record_list.append({
            'pattern': var_name,
            'var_name': var_name,
            # 'vaccine': rx_list[0]['vaccine_name'],
            # 'vaccine_type': rx_list[0]['vaccine_type'],
            'median_fold': float(median_fold),
            'num_ref_name': len(set([
                r['ref_name']
                for r in rx_list
                ])),
            'num_fold': num_fold,
            'S': num_s,
            'I': num_i,
            'R': num_r,
        })

        if var_name == 'other combo mut':
            rx_list += single_mut_records
            num_s, num_i, num_r, median_fold, num_fold = get_fold_stat(rx_list)
            record_list.append({
                'pattern': 'Individual mutations and mutation combinations',
                'var_name': 'Individual mutations and mutation combinations',
                # 'vaccine': rx_list[0]['vaccine_name'],
                # 'vaccine_type': rx_list[0]['vaccine_type'],
                'median_fold': '',
                'num_ref_name': len(set([
                    r['ref_name']
                    for r in rx_list
                    ])),
                'num_fold': num_fold,
                'S': '',
                'I': '',
                'R': '',
            })

    record_list = [
        i for i in record_list
        if i['var_name'] != 'other combo mut'
    ]

    record_list.sort(key=itemgetter(
        'var_name',
        'pattern',
        # 'vaccine',
        ))

    record_list = [
        i for i in record_list
        if not i['var_name'].startswith('Individual')
    ] + [
        i for i in record_list
        if i['var_name'].startswith('Individual')
    ]
    # record_list.append({
    #     'pattern': 'summary',
    #     'num_ref_name': len(set([
    #         r['ref_name']
    #         for r in db_records
    #     ])),
    #     'num_fold': sum([r['num_fold'] for r in db_records] + [0]),
    # })

    dump_csv(save_path, record_list)

    save_path = save_path.parent / (save_path.stem + '.json')
    dump_json(save_path, record_list)
