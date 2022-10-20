from importlib import import_module
import inspect
from preset import g


def is_task(func):
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    wrapper.is_task = True
    return wrapper


def check_func_is_task(func):
    if 'is_task' not in dir(func):
        return False

    return getattr(func, 'is_task')


def get_module_and_func(mdule_path):
    pass


def get_main_func(module_path):
    module_path_list = module_path.split('.')
    try:
        module = import_module(module_path)
        main_func_name = module_path_list[-1]
    except ModuleNotFoundError:
        module_path = '.'.join(module_path_list[:-1])
        main_func_name = module_path_list[-1]
        module = import_module(module_path)

    main_func = getattr(module, main_func_name)

    return main_func


def get_task_name(module_name):
    return module_name.split('.')[-1]


def get_task_info(task):
    if type(task) == str:
        func_name = task
        positional = []
        optional = {}
    elif type(task) == dict:
        func_name = list(task.keys())[0]
        arguments = list(task.values())[0]
        positional = [
            p
            for p in arguments
            if type(p) != dict
        ]
        optional = {
            k: v
            for p in arguments
            if type(p) == dict
            for k, v in p.items()
        }

    return (func_name, positional, optional)


def process_task(task):
    if type(task) == list:
        print('Work on task group', task)
        for sub_task in task:
            process_task(sub_task)
    else:
        func_name, positional, optional = get_task_info(task)
        print('Work on...', func_name, positional, optional)
        func = get_main_func(func_name)

        if test_require_conn(func):
            positional.insert(0, g.conn)

        func(*positional, **optional)


def test_require_conn(func):

    args = inspect.getargspec(func).args

    if args and 'conn' == args[0]:
        return True
    else:
        return False
