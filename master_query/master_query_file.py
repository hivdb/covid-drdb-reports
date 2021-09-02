import sqlite3
from pathlib import Path
from preset import load_csv, dump_csv

DB_PATH = Path('../covid-drdb/local/covid-drdb-latest.db')
DOMAINS = {
    'NTD': (1, 305),
    'RBD': (306, 534),
    'CTD': (535, 686),
    'S2': (687, 1273),
}

SPIKE_REF_SQL = """
SELECT
    position, amino_acid as refAA
FROM
    ref_amino_acid
WHERE
    gene = 'S'
ORDER BY
    position
;
"""


def run_query(sql_stat):
    conn = sqlite3.connect(str(DB_PATH))
    conn.row_factory = sqlite3.Row

    cursor = conn.cursor()
    cursor.execute(sql_stat)

    result = []
    for row in cursor.fetchall():
        rec = {}
        for key in row.keys():
            rec[key] = row[key]
        result.append(rec)
    return result

l = [i['position'] for i in run_query(SPIKE_REF_SQL)]

m = {i['position']: i['refAA'] for i in run_query(SPIKE_REF_SQL)}


DMS_SQL = """
SELECT
    *
FROM

"""

"""
select * from (SELECT rx_name, ab_name, count(ab_name) num_ab_name FROM 'rx_dms' group by rx_name) where num_ab_name = 1
"""



def query_dms_by_position(pos, aa):
    print('Position:', pos, ",AA:", aa)
    DMS_BINDING_SQL = """
    SELECT ace2_binding, expression FROM 'dms_ace2_binding' where position = {pos} and amino_acid = '{aa}'
    """.format(pos=pos, aa=aa)

    for i in run_query(DMS_BINDING_SQL):
        print(i)

    DMS_SCORE_SQL = """
    select a.ab_name, b.escape_score from (SELECT a.rx_name, b.ab_name FROM (select * from (SELECT rx_name, ab_name, count(ab_name) num_ab_name FROM 'rx_dms' group by rx_name) where num_ab_name = 1) a, 'antibodies' b where availability is not null and a.ab_name = b.ab_name) a, dms_escape_results b on a.rx_name = b.rx_name and  position = '{pos}' and amino_acid = '{aa}'
    """.format(pos=pos, aa=aa)

    result = run_query(DMS_SCORE_SQL)
    for i in result:
        show_values = []
        for v in i.values():
            show_values.append(str(v))
        print(':'.join(show_values))

def query_invitro_by_mutation(pos, aa):
    """
    SELECT b.ab_name, group_concat(a.ref_name, ',') FROM 'invitro_selection_results' a, rx_antibodies b on a.ref_name = b.ref_name and a.rx_name = b.rx_name and b.ab_name in ("ADG20","BRII-196","BRII-198","Bamlanivimab","C135","C144","Casirivimab","Cilgavimab","Etesevimab","Imdevimab","JMB2002","Regdanvimab","STI-2020","Sotrovimab","Tixagevimab","Vir-7832") and a.position = 484 and a.amino_acid = 'K' group by b.ab_name
    """
    pass

# def query_muation(pos, aa):
query_dms_by_position(445, 'A')