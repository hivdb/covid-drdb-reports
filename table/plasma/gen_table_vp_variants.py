from preset import DATA_FILE_PATH
from operator import itemgetter
from preset import dump_csv

from .common import gen_plasma_indiv_table
from .common import gen_plasma_aggre_table
from .common import convert_to_json
from variant.preset import KEY_VARIANTS
from .preset import INDIVIDUAL_VP_SQL
from .preset import AGGREGATED_VP_SQL
from .common import record_modifier


def gen_table_vp_variants(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_vp_variants.csv',
        json_save_path=DATA_FILE_PATH / 'table_vp_variants.json',
        ):
    indiv_records = gen_plasma_indiv_table(
        conn,
        KEY_VARIANTS,
        'rx_vacc_plasma',
        INDIVIDUAL_VP_SQL,
        plasma_type='VP',
        record_modifier=record_modifier
    )

    aggre_records = gen_plasma_aggre_table(
        conn,
        KEY_VARIANTS,
        'rx_vacc_plasma',
        AGGREGATED_VP_SQL,
        plasma_type='VP',
        record_modifier=record_modifier
    )

    records = indiv_records + aggre_records

    records.sort(key=itemgetter(
        'pattern', 'Plasma', 'ref_name'))

    dump_csv(csv_save_path, records)

    convert_to_json(json_save_path, records)
