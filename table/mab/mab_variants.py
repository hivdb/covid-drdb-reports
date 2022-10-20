from preset import DATA_FILE_PATH
from preset import dump_csv
from variant.preset import KEY_VARIANTS
from variant.preset import OMICRON_VARIANTS
from operator import itemgetter
from resistancy import RESISTANCE_FILTER
from resistancy import round_fold
from sql import get_sql_stmt_from_tmpl
from pathlib import Path
from preset import g

FOLDER = Path(__file__).parent.resolve()


def mab_variants(save_folder=DATA_FILE_PATH / 'mab'):

    _mab_variants(
        save_folder=save_folder,
        file_name='table_mab_variant',
        sql_tmpl=FOLDER / 'mab_variants.tmpl.sql',
        variants=KEY_VARIANTS,
        susc_table='susc_results_50_wt_view')

    _mab_variants(
        save_folder=save_folder,
        file_name='table_mab_omicron_variant',
        sql_tmpl=FOLDER / 'mab_variants.tmpl.sql',
        variants=OMICRON_VARIANTS,
        susc_table='susc_results_50_ba2_view')


def _mab_variants(save_folder, file_name, variants, sql_tmpl, susc_table):

    cursor = g.conn.cursor()

    records = []

    # uniq_record_list = set()

    for row_name, attr_r in variants.items():
        for _, resist_filter in RESISTANCE_FILTER.items():
            r_filter = attr_r.get('filter', [])
            filter = '\n    '.join(r_filter + resist_filter)

            sql = get_sql_stmt_from_tmpl(
                sql_tmpl, filters=filter, susc_table=susc_table)

            # print(sql)

            cursor.execute(sql)
            for row in cursor.fetchall():
                ref_name = row['ref_name']
                ab_name = row['ab_name']
                ab_class = row['class']

                fold = row['fold']
                # ineffective = row['ineffective']
                # if ineffective:
                #     fold = 1000
                fold = '{}'.format(round_fold(fold))

                iso_name = row_name
                if 'full genome' in iso_name:
                    ref_name = '{}*'.format(ref_name)
                    iso_name = iso_name.split()[0]

                # uniq_record = (
                #     iso_name,
                #     ab_name,
                #     ref_name,
                #     fold,
                # )
                # if uniq_record in uniq_record_list:
                #     continue
                # else:
                #     uniq_record_list.add(uniq_record)

                records.append({
                    'pattern': iso_name,
                    'ab_name': ab_name,
                    'class': ab_class or '',
                    # 'Resistance level': resist_name,
                    'fold': fold,
                    'ref_name': ref_name
                })

    records.sort(key=itemgetter(
        'pattern',
        'class',
        'ab_name'))

    dump_csv(save_folder / f'{file_name}.csv', records)

    # json_records = defaultdict(list)
    # for r in records:
    #     variant = r['pattern']
    #     json_records[variant].append({
    #         'variant': variant,
    #         'rx_name': r['ab_name'],
    #         'mab_class': r['class'],
    #         # 'fold': r['Fold'].replace('>', '&gt;'),
    #         'fold': r['fold'],
    #         'ref_name': r['ref_name']
    #     })

    # records = []
    # for variant, assays in json_records.items():
    #     records.append({
    #         'variant': variant,
    #         'assays': sorted(assays, key=itemgetter('mab_class')),
    #     })

    # variant = sorted(records, key=itemgetter('variant'))
