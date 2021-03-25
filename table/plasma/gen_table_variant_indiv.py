from preset import DATA_FILE_PATH

from .common import gen_plasma_table
from .prefix import SUBROWS
from .prefix import VARIANTS
from .prefix import INDIVIDUAL_RESULTS_SQL
from .common import record_modifier


def gen_table_plasma_variant(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_plasma_variant.csv',
        json_save_path=DATA_FILE_PATH / 'table_plasma_variant.json',
        ):
    gen_plasma_table(
        conn,
        VARIANTS,
        SUBROWS,
        csv_save_path,
        json_save_path,
        INDIVIDUAL_RESULTS_SQL,
        record_modifier
    )
