from preset import DATA_FILE_PATH
from preset import dump_csv
from variant.preset import ONE_MUT_VARIANT
from variant.preset import COMBO_MUT_VARIANT
from variant.preset import group_by_variant
from plasma.preset import AGGREGATED_RESULTS_SQL
from variant.preset import CONTROL_VARIANTS_SQL


RXTYPE = {
    'rx_vacc_plasma': DATA_FILE_PATH / 'summary_vp_aggre.csv',
    'rx_conv_plasma': DATA_FILE_PATH / 'summary_cp_aggre.csv',
}


def gen_table_variant_aggre(conn):

    cursor = conn.cursor()
    for rxtype, save_path in RXTYPE.items():
        sql = AGGREGATED_RESULTS_SQL.format(
            rxtype=rxtype,
            filters='',
            control_variants=CONTROL_VARIANTS_SQL,
            )

        # print(sql)
        cursor.execute(sql)

        records = cursor.fetchall()

        indiv_records, combo_records = group_by_variant(records)

        results = []
        results += get_fold_results(indiv_records, 'indiv')
        results += get_fold_results(combo_records, 'combo')

        dump_csv(save_path, results)


def get_variant_group(group_name):
    if group_name == 'indiv':
        return ONE_MUT_VARIANT
    else:
        return COMBO_MUT_VARIANT


def get_fold_results(mut_group, variant_group_name):
    results = []
    variant_group = get_variant_group(variant_group_name)
    for variant, r_list in mut_group.items():
        for r in r_list:
            variant_info = variant_group.get(variant)
            results.append({
                'pattern': variant,
                'domain': variant_info.get('domain'),
                'varname': variant_info.get('varname'),
                'reference': r['ref_name'],
                'fold_cmp': r['fold_cmp'],
                'median': r['fold'],
                'count': r['sample_count']
            })

    return results
