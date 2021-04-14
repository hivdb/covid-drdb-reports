from preset import DATA_FILE_PATH
from operator import itemgetter
from preset import dump_csv

from .common import gen_plasma_indiv_table
from .common import gen_plasma_aggre_table
from .common import convert_to_json
from .preset import SUBROWS
from .preset import MUTATIONS
from .preset import INDIVIDUAL_RESULTS_SQL
from .preset import AGGREGATED_RESULTS_SQL


def gen_table_plasma_muts(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_plasma_muts.csv',
        json_save_path=DATA_FILE_PATH / 'table_plasma_muts.json'
        ):
    indiv_records = gen_plasma_indiv_table(
        conn,
        MUTATIONS,
        SUBROWS,
        INDIVIDUAL_RESULTS_SQL
    )

    aggre_records = gen_plasma_aggre_table(
        conn,
        MUTATIONS,
        SUBROWS,
        AGGREGATED_RESULTS_SQL
    )

    records = indiv_records + aggre_records

    records.sort(key=itemgetter(
        'Variant name', 'Plasma', 'Reference'))

    dump_csv(csv_save_path, records)

    convert_to_json(json_save_path, records)
