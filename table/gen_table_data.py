import sys
import sqlite3
from pathlib import Path

from variant.gen_table_key_variant import gen_table_key_variant
from variant.preset import get_grouped_variants
from variant.gen_table_variant import gen_table_variant

from mab.preset import init_synonyms_map
from mab.preset import init_abname2class
from mab.gen_table_mab import gen_table_mab
from mab.gen_table_mab_variant import gen_table_mab_variant
from mab.gen_table_mab_muts import gen_table_mab_muts
from mab.gen_table_all_mab import gen_table_all_mab

from plasma.gen_table_summary import gen_table_plasma
from plasma.gen_table_variant_indiv import gen_table_plasma_variant
from plasma.gen_table_mut_indiv import gen_table_plasma_muts
from plasma.gen_table_variant_aggre import gen_table_plasma_variant_aggre
from plasma.gen_table_mut_aggre import gen_table_plasma_muts_aggre
from plasma.gen_table_vp import gen_table_vp

from lookup_view import get_aggregated_studies


def gen_report(db_path):
    db_path = Path(db_path).resolve()
    print(db_path)
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row

    init_synonyms_map(conn)
    init_abname2class(conn)

    get_grouped_variants(conn)

    get_aggregated_studies(conn)

    gen_table_plasma_variant(conn)
    gen_table_plasma_variant_aggre(conn)
    gen_table_plasma_muts(conn)
    gen_table_plasma_muts_aggre(conn)
    gen_table_vp(conn)

    gen_table_mab_variant(conn)
    gen_table_mab_muts(conn)
    gen_table_all_mab(conn)

    gen_table_key_variant(conn)
    gen_table_variant(conn)

    gen_table_plasma()
    gen_table_mab()

    print('done')


if __name__ == '__main__':
    db_path = sys.argv[1]
    gen_report(db_path)
