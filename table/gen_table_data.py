import sys
import sqlite3
from pathlib import Path
from preset import init_synonyms_map
from preset import init_abname2class


from gen_table_summary import gen_table_summary
from gen_table_mab import gen_table_mab
from gen_table_mab_variant import gen_table_mab_variant
from gen_table_mab_muts import gen_table_mab_muts

from plasma.gen_table_summary import gen_table_plasma
from plasma.gen_table_variant_indiv import gen_table_plasma_variant
from plasma.gen_table_mut_indiv import gen_table_plasma_muts
from plasma.gen_table_variant_aggre import gen_table_plasma_variant_aggre
from plasma.gen_table_mut_aggre import gen_table_plasma_muts_aggre
from plasma.gen_table_vp import gen_table_vp

from gen_escape import gen_escape
from gen_binding import gen_table_binding
from gen_dms import gen_dms

from lookup_view import get_aggregated_studies


def gen_report(db_path):
    db_path = Path(db_path).resolve()
    print(db_path)
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    init_synonyms_map(conn)
    init_abname2class(conn)

    get_aggregated_studies(conn)

    gen_table_plasma_variant(conn)
    gen_table_plasma_variant_aggre(conn)
    gen_table_plasma_muts(conn)
    gen_table_plasma_muts_aggre(conn)
    gen_table_vp(conn)

    gen_table_mab_variant(conn)
    gen_table_mab_muts(conn)

    gen_table_summary(conn)
    gen_table_plasma()
    gen_table_mab()

    gen_escape(conn)
    gen_table_binding(conn)
    gen_dms(conn)
    print('done')


if __name__ == '__main__':
    db_path = sys.argv[1]
    gen_report(db_path)
