#!/usr/bin/env python3
from preset import DATA_FILE_PATH
from preset import init_db
from task import process_task
import click


@click.command()
@click.argument('db_path')
def work(db_path, save_path=DATA_FILE_PATH / 'ref_names.md'):

    conn = init_db(db_path)

    cursor = conn.cursor()

    cursor.execute(
        'SELECT distinct ref_name from treatments order by ref_name;')

    with open(save_path, 'w') as fd:
        for i in cursor.fetchall():
            fd.write(f"1. **{i['ref_name']}**: [^{i['ref_name']}#inline]\n")


if __name__ == '__main__':
    print('working...')
    work()
    print('done.')
