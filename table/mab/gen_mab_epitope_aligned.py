from preset import DATA_FILE_PATH
from preset import dump_csv
from sql import row2dict
from .preset import ACE2_FOOT_PRINT_WT


SUMMARY_SQL = """
SELECT DISTINCT
    mab.*,
    epi.epitope
FROM
    single_mab_view mab,
    mab_epitope_view epi
WHERE
    mab.ab_name = epi.ab_name
    AND
    mab.availability IS NOT NULL
;
"""

REF_AA_SQL = """
SELECT
    position,
    amino_acid
FROM
    ref_amino_acid
WHERE
    gene = 'S'
    AND
    position >= 331
    AND
    position <= 550
;
"""


def gen_mab_epitope_aligned(
        conn,
        save_path=DATA_FILE_PATH / 'mab' / 'main_mab_epitope_aligned.csv'):

    cursor = conn.cursor()

    sql = SUMMARY_SQL

    cursor.execute(sql)

    records = row2dict(cursor.fetchall())
    dump_csv(DATA_FILE_PATH / 'mab' / 'main_mab_epitope.csv', records)

    pos_ref = get_pos_ref(conn)
    headers = sorted(pos_ref.keys())
    headers = ['mab'] + headers
    first_row = {
        'mab': ''
    }
    for pos in sorted(pos_ref.keys()):
        first_row[pos] = pos_ref[pos]

    second_row = {
        'mab': 'ACE2_footprint',
    }
    for pos in sorted(pos_ref.keys()):
        if pos not in ACE2_FOOT_PRINT_WT:
            continue
        second_row[pos] = '*'

    results = []
    results.append(first_row)
    results.append(second_row)
    for rec in records:
        mab = rec['ab_name']
        epitope = rec['epitope']
        epi_list = [int(i) for i in epitope.split('+')]
        new_rec = {
            'mab': mab
        }
        for pos in epi_list:
            new_rec[pos] = '*'

        results.append(new_rec)

    dump_csv(save_path, results)


def get_pos_ref(conn):
    cursor = conn.cursor()

    sql = REF_AA_SQL

    cursor.execute(sql)

    results = row2dict(cursor.fetchall())

    pos_ref = {}
    for row in results:
        ref_aa = row['amino_acid']
        pos_aa = row['position']
        pos_ref[pos_aa] = ref_aa

    return pos_ref
