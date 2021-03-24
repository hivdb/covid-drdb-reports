import sys
import sqlite3
from preset import init_synonyms_map
from preset import init_abname2class


from gen_table_summary import gen_table_summary
from gen_table_plasma import gen_table_plasma
from gen_table_mab import gen_table_mab
from gen_table_plasma_variant import gen_table_plasma_variant
from gen_table_mab_variant import gen_table_mab_variant
from gen_table_plasma_muts import gen_table_plasma_muts
from gen_table_mab_muts import gen_table_mab_muts


def gen_report(db_path):
    conn = sqlite3.connect(db_path)
    init_synonyms_map(conn)
    init_abname2class(conn)

    gen_table_plasma_variant(conn)
    gen_table_mab_variant(conn)
    gen_table_plasma_muts(conn)
    gen_table_mab_muts(conn)

    gen_table_summary(conn)
    gen_table_plasma()
    gen_table_mab()
    print('done')


if __name__ == '__main__':
    db_path = sys.argv[1]
    gen_report(db_path)
