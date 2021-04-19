from collections import defaultdict
from statistics import median
from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter
from .preset import INDIV_VARIANT
from .preset import COMBO_VARIANT
from preset import round_number
from resistancy import is_partial_resistant
from resistancy import is_resistant
from resistancy import is_susc
from variant.preset import CONTROL_VARIANTS_SQL

SQL = """
SELECT
    r.vaccine_name,
    s.variant_name,
    s.fold,
    s.cumulative_count as count
FROM
    susc_results as s
INNER JOIN rx_immu_plasma as r ON
    s.ref_name = r.ref_name
    AND s.rx_name = r.rx_name
WHERE s.control_variant_name in {control_variants}
    AND s.fold IS NOT NULL;
""".format(control_variants=CONTROL_VARIANTS_SQL)


def gen_table_variant_vp(conn):
    by_variant(
        conn,
        'indiv',
        save_path=DATA_FILE_PATH / 'summary_variant_indiv_vp.csv'
    )
    by_variant(
        conn,
        'combo',
        save_path=DATA_FILE_PATH / 'summary_variant_combo_vp.csv'
    )


def by_variant(conn, indiv_or_combo, save_path):
    if indiv_or_combo == 'indiv':
        variant_mapper = INDIV_VARIANT
    else:
        variant_mapper = COMBO_VARIANT
    cursor = conn.cursor()

    cursor.execute(SQL)

    variant_group = defaultdict(list)
    for rec in cursor.fetchall():
        variant = rec['variant_name']
        variant = variant_mapper.get(variant)
        if not variant:
            continue
        variant = variant['disp']
        variant_group[variant].append(rec)

    record_list = []
    for variant, rlist in variant_group.items():
        vacc_group = defaultdict(list)
        for rec in rlist:
            vacc_name = rec['vaccine_name']
            vacc_group[vacc_name].append(rec)

        for vacc_name, rx_list in vacc_group.items():
            all_fold = [(r['fold'], r['count']) for r in rx_list if r['fold']]
            num_s = sum([r[1] for r in all_fold if is_susc(r[0])])
            num_i = sum([r[1] for r in all_fold if is_partial_resistant(r[0])])
            num_r = sum([r[1] for r in all_fold if is_resistant(r[0])])

            all_fold = [[i[0]] * i[1] for i in all_fold]
            all_fold = [i for j in all_fold for i in j if i]
            median_fold = round_number(median(all_fold)) if all_fold else ''

            num_results = sum([r['count'] for r in rx_list] + [0])

            if indiv_or_combo == 'indiv':
                variant_info = INDIV_VARIANT.get(variant)
                record_list.append({
                    'variant': variant,
                    'vaccine': vacc_name,
                    'RefAA': variant_info['ref_aa'],
                    'Position': variant_info['position'],
                    'AA': variant_info['aa'],
                    'Domain': variant_info['domain'],
                    'median_fold': median_fold,
                    'samples': num_results,
                    'S': num_s,
                    'I': num_i,
                    'R': num_r,
                })
            else:
                variant_info = COMBO_VARIANT.get(variant)
                nickname = variant_info['nickname']
                record_list.append({
                    'variant': variant,
                    'nickname': nickname,
                    'vaccine': vacc_name,
                    'median_fold': median_fold,
                    'samples': num_results,
                    'S': num_s,
                    'I': num_i,
                    'R': num_r,
                })

    record_list.sort(key=itemgetter(
        'variant',
        'vaccine',
        ))

    dump_csv(save_path, record_list)
