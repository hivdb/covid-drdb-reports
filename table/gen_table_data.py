import sys
import sqlite3
from preset import init_synonyms_map
from preset import init_abname2class
from gen_table3 import gen_table3
from gen_tableS1 import gen_tableS1
from gen_tableS2 import gen_tableS2
from gen_tableS3 import gen_tableS3
from gen_tableS4 import gen_tableS4


def gen_report(db_path):
    conn = sqlite3.connect(db_path)
    init_synonyms_map(conn)
    init_abname2class(conn)

    gen_table3(conn)
    gen_tableS1(conn)
    gen_tableS2(conn)
    gen_tableS3(conn)
    gen_tableS4(conn)


if __name__ == '__main__':
    db_path = sys.argv[1]
    gen_report(db_path)

