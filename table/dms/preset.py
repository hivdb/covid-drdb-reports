

DMS_POSITIONS = []
DMS_POS_SQL = """
SELECT
    distinct position
FROM
    dms_escape_results;
"""


def get_dms_positions(conn):
    global DMS_POSITIONS

    cursor = conn.cursor()
    cursor.execute(DMS_POS_SQL)

    records = cursor.fetchall()
    DMS_POSITIONS = [rec['position'] for rec in records]
