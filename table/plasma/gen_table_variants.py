from preset import DATA_FILE_PATH
from operator import itemgetter
from preset import dump_csv

from .common import gen_plasma_indiv_table
from .common import gen_plasma_aggre_table
from .common import convert_to_json
from .preset import SUBROWS
from .preset import VARIANTS
from .preset import INDIVIDUAL_RESULTS_SQL
from .preset import AGGREGATED_RESULTS_SQL
from .common import record_modifier


def gen_table_plasma_variant(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_plasma_variant.csv',
        json_save_path=DATA_FILE_PATH / 'table_plasma_variant.json',
        ):
    indiv_records = gen_plasma_indiv_table(
        conn,
        VARIANTS,
        SUBROWS,
        INDIVIDUAL_RESULTS_SQL,
        record_modifier
    )

    aggre_records = gen_plasma_aggre_table(
        conn,
        VARIANTS,
        SUBROWS,
        AGGREGATED_RESULTS_SQL,
        record_modifier
    )

    records = indiv_records + aggre_records

    records.sort(key=itemgetter(
        'Variant name', 'Plasma', 'Reference'))

    dump_csv(csv_save_path, records)

    convert_to_json(json_save_path, records)
