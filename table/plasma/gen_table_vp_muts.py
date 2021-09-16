from preset import DATA_FILE_PATH
from operator import itemgetter
from preset import dump_csv

from .common import gen_plasma_indiv_table
from .common import gen_plasma_aggre_table
from .common import convert_to_json
from .preset import VP_FILTER
from .preset import MUTATIONS
from .preset import INDIVIDUAL_VP_SQL
from .preset import AGGREGATED_VP_SQL


def gen_table_vp_muts(
        conn,
        csv_save_path=DATA_FILE_PATH / 'table_vp_muts.csv',
        json_save_path=DATA_FILE_PATH / 'table_vp_muts.json'
        ):
    indiv_records = gen_plasma_indiv_table(
        conn,
        MUTATIONS,
        VP_FILTER,
        INDIVIDUAL_VP_SQL,
        plasma_type='VP',
    )

    aggre_records = gen_plasma_aggre_table(
        conn,
        MUTATIONS,
        VP_FILTER,
        AGGREGATED_VP_SQL,
        plasma_type='VP',
    )

    records = indiv_records + aggre_records

    records.sort(key=itemgetter(
        'pattern', 'Plasma', 'Reference'))

    dump_csv(csv_save_path, records)

    convert_to_json(json_save_path, records)
