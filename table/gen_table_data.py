import sys
import sqlite3
from pathlib import Path

from variant.gen_table_key_variant import gen_table_key_variant
from variant.preset import get_grouped_variants
from variant.preset import get_spike_ref
from variant.gen_table_variant import gen_table_variant
from variant.gen_table_variant_mab import gen_table_variant_mab
from variant.gen_table_variant_vp import gen_table_variant_vp
from variant.gen_table_variant_cp import gen_table_variant_cp
from variant.gen_table_variant_aggre import gen_table_variant_aggre

from mab.gen_table_mab import gen_table_mab
from mab.gen_table_mab_variant import gen_table_mab_variant
from mab.gen_table_mab_muts import gen_table_mab_muts
from mab.gen_table_all_mab import gen_table_all_mab

from plasma.gen_table_summary import gen_table_plasma
from plasma.gen_table_variants import gen_table_plasma_variant
from plasma.gen_table_muts import gen_table_plasma_muts
from plasma.gen_table_vp import gen_table_vp

from study.gen_study import gen_study


def gen_report(db_path):
    db_path = Path(db_path).resolve()
    print(db_path)
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row

    get_spike_ref(conn)

    get_grouped_variants(conn)

    gen_table_plasma_variant(conn)
    gen_table_plasma_muts(conn)
    gen_table_vp(conn)

    gen_table_mab_variant(conn)
    gen_table_mab_muts(conn)
    gen_table_all_mab(conn)

    gen_table_key_variant(conn)
    gen_table_variant(conn)
    gen_table_variant_mab(conn)
    gen_table_variant_vp(conn)
    gen_table_variant_cp(conn)
    gen_table_variant_aggre(conn)

    gen_study(conn)

    gen_table_plasma()
    gen_table_mab()

    print('done')


if __name__ == '__main__':
    db_path = sys.argv[1]
    gen_report(db_path)
