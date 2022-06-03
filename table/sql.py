def get_sql_stmt(file_path):
    with open(file_path) as fd:
        stmt = fd.read()

    stmt = stmt.strip().rstrip(';')
    return stmt


def get_sql_stmt_from_tmpl(file_path, **kargs):
    with open(file_path) as fd:
        stmt_tmpl = fd.read()

    stmt = stmt_tmpl.format(**kargs)
    stmt = stmt.strip().rstrip(';')
    return stmt


def row2dict(rows):
    result = []
    for row in rows:
        rec = {}
        for key in row.keys():
            rec[key] = row[key]
        result.append(rec)

    return result
