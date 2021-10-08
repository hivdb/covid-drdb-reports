from preset import DATA_FILE_PATH
from operator import itemgetter
from preset import dump_csv

from .common import gen_plasma_indiv_table
from .common import gen_plasma_aggre_table
from .common import convert_to_json
from variant.preset import KEY_MUTATIONS
from .preset import INDIVIDUAL_CP_SQL
from .preset import AGGREGATED_CP_SQL


def gen_table_cp_muts(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_cp_muts.csv',
        json_save_path=DATA_FILE_PATH / 'table_cp_muts.json'
        ):
    indiv_records = gen_plasma_indiv_table(
        conn,
        KEY_MUTATIONS,
        'rx_conv_plasma_wt_infect_view',
        INDIVIDUAL_CP_SQL,
        'CP',
    )

    aggre_records = gen_plasma_aggre_table(
        conn,
        KEY_MUTATIONS,
        'rx_conv_plasma_wt_infect_view',
        AGGREGATED_CP_SQL,
        'CP',
    )

    records = indiv_records + aggre_records

    records.sort(key=itemgetter(
        'pattern', 'Plasma', 'Reference'))

    dump_csv(csv_save_path, records)

    convert_to_json(json_save_path, records)
