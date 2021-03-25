from preset import DATA_FILE_PATH

from .common import gen_plasma_aggre_table
from .prefix import SUBROWS
from .prefix import MUTATIONS
from .prefix import AGGREGATED_RESULTS_SQL


def gen_table_plasma_muts_aggre(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_plasma_muts_aggre.csv',
        json_save_path=DATA_FILE_PATH / 'table_plasma_muts_aggre.json'
        ):
    gen_plasma_aggre_table(
        conn,
        MUTATIONS,
        SUBROWS,
        csv_save_path,
        json_save_path,
        AGGREGATED_RESULTS_SQL
    )
