from preset import DATA_FILE_PATH

from .common import gen_plasma_aggre_table
from .prefix import SUBROWS
from .prefix import VARIANTS
from .prefix import AGGREGATED_RESULTS_SQL
from .common import record_modifier


def gen_table_plasma_variant_aggre(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_plasma_variant_aggre.csv',
        json_save_path=DATA_FILE_PATH / 'table_plasma_variant_aggre.json',
        ):
    gen_plasma_aggre_table(
        conn,
        VARIANTS,
        SUBROWS,
        csv_save_path,
        json_save_path,
        AGGREGATED_RESULTS_SQL,
        record_modifier
    )
