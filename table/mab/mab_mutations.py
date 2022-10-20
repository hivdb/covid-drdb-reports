from preset import DATA_FILE_PATH
from preset import dump_csv
from pathlib import Path
from variant.preset import KEY_MUTATIONS
from variant.preset import OMICRON_MUTATIONS
from operator import itemgetter
from collections import defaultdict
from resistancy import RESISTANCE_FILTER
from resistancy import round_fold
from preset import g
from sql import get_sql_stmt_from_tmpl


FOLDER = Path(__file__).parent.resolve()


def mab_mutations(save_folder=DATA_FILE_PATH / 'mab'):

    _mab_mutations(
        save_folder=save_folder,
        file_name='table_mab_muts',
        sql_tmpl=FOLDER / 'mab_mutations.tmpl.sql',
        mutations=KEY_MUTATIONS,
        susc_table='susc_results_50_wt_view')

    _mab_mutations(
        save_folder=save_folder,
        file_name='table_mab_omicron_muts',
        sql_tmpl=FOLDER / 'mab_mutations.tmpl.sql',
        mutations=OMICRON_MUTATIONS,
        susc_table='susc_results_50_ba2_view')


def _mab_mutations(
        save_folder,
        file_name,
        sql_tmpl,
        mutations,
        susc_table):
    cursor = g.conn.cursor()

    records = []
    for row_name, attr_r in mutations.items():
        for _, resist_filter in RESISTANCE_FILTER.items():
            iso_type = attr_r.get('iso_type')

            r_filter = attr_r.get('filter', [])
            filter = '\n    '.join(r_filter + resist_filter)

            sql = get_sql_stmt_from_tmpl(
                sql_tmpl,
                iso_type=iso_type,
                filters=filter,
                susc_table=susc_table)

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

                records.append({
                    'pattern': row_name,
                    'ab_name': ab_name,
                    'class': ab_class or '',
                    # 'Resistance level': resist_name,
                    'fold': fold,
                    'ref_name': ref_name
                })

    records.sort(key=itemgetter(
        'pattern',
        'class',
        'ab_name',
        ))

    dump_csv(save_folder / f'{file_name}.csv', records)

    # json_records = defaultdict(list)
    # for r in records:
    #     variant = r['pattern']
    #     json_records[variant].append({
    #         'variant': variant,
    #         'rx_name': r['ab_name'],
    #         'mab_class': r['class'],
    #         # 'fold': r['fold'].replace('>', '&gt;'),
    #         'fold': r['fold'],
    #         'ref_name': r['ref_name']
    #     })

    # records = []
    # for pattern, assays in json_records.items():
    #     records.append({
    #         'pattern': pattern,
    #         'assays': sorted(assays, key=itemgetter('mab_class')),
    #     })

    # records.sort(key=itemgetter('pattern'))
