from preset import DATA_FILE_PATH
from operator import itemgetter
from preset import dump_csv

from .common import gen_plasma_indiv_table
from .common import gen_plasma_aggre_table
from .common import convert_to_json
from .preset import VP_FILTER
from .preset import VARIANTS
from .preset import INDIVIDUAL_RESULTS_SQL
from .preset import AGGREGATED_RESULTS_SQL
from .common import record_modifier


def gen_table_vp_variants(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_vp_variants.csv',
        json_save_path=DATA_FILE_PATH / 'table_vp_variants.json',
        ):
    indiv_records = gen_plasma_indiv_table(
        conn,
        VARIANTS,
        VP_FILTER,
        INDIVIDUAL_RESULTS_SQL,
        record_modifier=record_modifier
    )

    aggre_records = gen_plasma_aggre_table(
        conn,
        VARIANTS,
        VP_FILTER,
        AGGREGATED_RESULTS_SQL,
        record_modifier=record_modifier
    )

    records = indiv_records + aggre_records

    records.sort(key=itemgetter(
        'Variant name', 'Plasma', 'Reference'))

    dump_csv(csv_save_path, records)

    convert_to_json(json_save_path, records)
