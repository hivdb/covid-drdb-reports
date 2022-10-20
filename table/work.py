#!/usr/bin/env python3
from pathlib import Path
from preset import load_yaml
from preset import init_db
from task import process_task
import click


@click.command()
@click.argument('db_path')
@click.argument('task_name', required=False)
def work(db_path, task_name, task_file='tasks.yml'):
    init_db(db_path)

    tasks = load_yaml(Path(__file__).resolve().parent / task_file)

    main_tasks = tasks['main_tasks']
    if task_name:
        main_tasks = [task_name]

    assert (type(main_tasks) == list), 'main task error'

    task_groups = tasks['task_groups']

    main_tasks = [
        task_groups[t] if (type(t) == str and t) in task_groups else t
        for t in main_tasks
    ]

    for task in main_tasks:
        process_task(task)


if __name__ == '__main__':
    print('working...')
    work()
    print('done.')
