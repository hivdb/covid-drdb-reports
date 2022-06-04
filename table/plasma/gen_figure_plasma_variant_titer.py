from sql import get_sql_stmt
from sql import get_sql_stmt_from_tmpl
from preset import DATA_FILE_PATH
from preset import dump_csv
from pathlib import Path
from sql import row2dict

FOLDER = Path(__file__).parent.resolve()


def gen_figure_plasma_variant_titer(
        conn,
        save_path=DATA_FILE_PATH / 'figure' / 'figure_plasma_variant_titer.csv'):

    vp_plasma_history = get_sql_stmt(FOLDER / 'vp_plasma_history.sql')
    plasma_summary = get_sql_stmt_from_tmpl(
        FOLDER / 'plasma_summary.sql', vp_plasma_history=vp_plasma_history)

    sql = get_sql_stmt_from_tmpl(
        FOLDER / 'plasma_variant_titer.tmpl.sql',
        plasma_summary=plasma_summary)

    table = conn.execute(sql)
    table = row2dict(table)

    [
        i.update({
            'month': max(int(i['days'] / 30), 1),
            'vaccine_name':
                ' + '.join(
                    sorted(list(set(
                        j for j in (i.get('vaccine_names') or '').split('+')
                        if j
                        )))
                ),
            'infection':
                ' + '.join(
                    sorted(list(set(
                        j for j in (i.get('infections') or '').split('+')
                        if j
                        )))
                )
        })
        for i in table
    ]

    dump_csv(save_path, table)
