from preset import DATA_FILE_PATH
from operator import itemgetter
from preset import dump_csv

from .common import gen_plasma_indiv_table
from .common import gen_plasma_aggre_table
from variant.preset import KEY_VARIANTS
from .preset import INDIVIDUAL_CP_SQL
from .preset import AGGREGATED_CP_SQL
from .common import record_modifier


def gen_table_cp_variants(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_cp_variants.csv',
        ):
    indiv_records = gen_plasma_indiv_table(
        conn,
        KEY_VARIANTS,
        'rx_conv_plasma_wt_infect_view',
        INDIVIDUAL_CP_SQL,
        'CP',
        record_modifier=record_modifier
    )

    aggre_records = gen_plasma_aggre_table(
        conn,
        KEY_VARIANTS,
        'rx_conv_plasma_wt_infect_view',
        AGGREGATED_CP_SQL,
        'CP',
        record_modifier=record_modifier
    )

    records = indiv_records + aggre_records

    records.sort(key=itemgetter(
        'pattern', 'Plasma', 'ref_name'))

    dump_csv(csv_save_path, records)
