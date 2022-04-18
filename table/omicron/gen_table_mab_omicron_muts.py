from preset import DATA_FILE_PATH
from preset import dump_csv
from variant.preset import OMICRON_MUTATIONS
from mab.preset import MAIN_MAB
from operator import itemgetter
from resistancy import RESISTANCE_FILTER
from resistancy import round_fold
from preset import group_records_by


MAB_MUTS_SQL = """
SELECT
    s.ref_name,
    rx.ab_name,
    rx.class,
    mut.position,
    mut.domain,
    s.fold_cmp,
    s.fold,
    -- control_pot.potency as control_ic50,
    -- test_pot.potency as test_ic50,
    s.ineffective
FROM
    susc_results_50_wt_view as s,
    rx_mab_view as rx,
    -- rx_potency as control_pot,
    -- rx_potency as test_pot,
    {iso_type} mut
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name

    -- AND
    -- s.ref_name = control_pot.ref_name
    -- AND
    -- s.rx_name = control_pot.rx_name
    -- AND
    -- s.control_iso_name = control_pot.iso_name

    -- AND
    -- s.ref_name = test_pot.ref_name
    -- AND
    -- s.rx_name = test_pot.rx_name
    -- AND
    -- s.iso_name = test_pot.iso_name

    AND
    s.iso_name = mut.iso_name
    AND
    s.fold IS NOT NULL
    AND (
        rx.availability IS NOT NULL
        OR rx.pdb_id IS NOT NULL
    )
    AND
    {filters}
"""


def gen_table_mab_omicron_muts(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_mab_omicron_muts.csv',
        ):
    cursor = conn.cursor()

    records = []
    for row_name, attr_r in OMICRON_MUTATIONS.items():
        for resist_name, resist_filter in RESISTANCE_FILTER.items():
            iso_type = attr_r.get('iso_type')

            r_filter = attr_r.get('filter', [])
            filter = '\n    '.join(r_filter + resist_filter)
            sql = MAB_MUTS_SQL.format(
                iso_type=iso_type,
                filters=filter,
            )
            # print(sql)

            cursor.execute(sql)
            for row in cursor.fetchall():
                ref_name = row['ref_name']
                ab_name = row['ab_name']
                ab_class = row['class']
                position = row['position']
                domain = row['domain']
                # control_ic50 = row['control_ic50']
                # test_ic50 = row['test_ic50']

                fold = row['fold']
                # ineffective = row['ineffective']
                # if ineffective:
                #     fold = 1000
                fold = '{}'.format(round_fold(fold))

                records.append({
                    'pattern': row_name,
                    'position': position,
                    'domain': domain,
                    'ab_name': ab_name,
                    'class': ab_class or '',
                    'fold': fold,
                    # 'control_ic50': control_ic50,
                    # 'test_ic50': test_ic50,
                    'ref_name': ref_name
                })

    records.sort(key=itemgetter(
        'position',
        'class',
        'ab_name',
        ))

    dump_csv(csv_save_path, records)

    get_main_mab_fold(records)
    # get_main_mab_ic50(records)


def get_main_mab_fold(
        records,
        save_path=DATA_FILE_PATH / 'omicron' / 'main_mab_omicron_muts_fold.csv'):

    mut_group = group_records_by(records, 'pattern')

    results = []
    for mut_pattern, rec_list in mut_group.items():
        position = rec_list[0]['position']
        domain = rec_list[0]['domain']

        rec_list = [i for i in rec_list if i['ab_name'] in MAIN_MAB.keys()]

        ref_group = group_records_by(rec_list, 'ref_name')

        for ref_name, ref_rec_list in ref_group.items():

            rec = {
                'aa_mut': mut_pattern,
                'position': position,
                'domain': domain,
                'ref_name': ref_name,
            }
            for mab_name, acronym in MAIN_MAB.items():
                _ref_rec_list = [
                    i for i in ref_rec_list if i['ab_name'] == mab_name]
                if not _ref_rec_list:
                    rec[acronym] = ''
                    continue
                for _rec in _ref_rec_list:
                    rec[acronym] = _rec['fold']

            results.append(rec)

    dump_csv(save_path, results)


def get_main_mab_ic50(
        records,
        save_path=DATA_FILE_PATH / 'omicron' / 'main_mab_omicron_muts_ic50.csv'):

    results = []
    for rec in records:
        del rec['class']
        if rec['ab_name'] not in MAIN_MAB.keys():
            continue

        rec['mab'] = MAIN_MAB[rec['ab_name']]
        del rec['ab_name']

        results.append(rec)

    dump_csv(save_path, results)
