from preset import DATA_FILE_PATH

from .common import gen_plasma_table
from .preset import SUBROWS
from .preset import MUTATIONS
from .preset import INDIVIDUAL_RESULTS_SQL


def gen_table_plasma_muts(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_plasma_muts.csv',
        json_save_path=DATA_FILE_PATH / 'table_plasma_muts.json'
        ):
    gen_plasma_table(
        conn,
        MUTATIONS,
        SUBROWS,
        csv_save_path,
        json_save_path,
        INDIVIDUAL_RESULTS_SQL
    )
