import sys
import sqlite3
from preset import init_synonyms_map
from preset import init_abname2class
from gen_table3 import gen_table3
from gen_table4 import gen_table4
from gen_tableS3 import gen_tableS3
from gen_tableS4 import gen_tableS4
from gen_tableS5 import gen_tableS5
from gen_tableS6 import gen_tableS6
from gen_tableS7 import gen_tableS7


def gen_report(db_path):
    conn = sqlite3.connect(db_path)
    init_synonyms_map(conn)
    init_abname2class(conn)

    gen_tableS3(conn)
    gen_tableS4(conn)
    gen_tableS5(conn)
    gen_tableS6(conn)
    gen_tableS7(conn)

    gen_table3()
    gen_table4()


if __name__ == '__main__':
    db_path = sys.argv[1]
    gen_report(db_path)
