import csv
from pathlib import Path
import decimal
from decimal import Decimal
import json
from decimal import localcontext


WS = Path(__file__).absolute().parent.parent
DATA_FILE_PATH = WS / 'report_tables'
DATA_FILE_PATH.mkdir(exist_ok=True)


def dump_csv(file_path, records, headers=[]):
    if not records:
        return
    if not headers and records:
        headers = records[0].keys()

    with open(file_path, 'w', encoding='utf-8-sig') as fd:
        writer = csv.DictWriter(fd, fieldnames=headers)
        writer.writeheader()
        writer.writerows(records)


def load_csv(file_path):
    records = []
    with open(file_path, encoding='utf-8-sig') as fd:
        for record in csv.DictReader(fd):
            records.append(record)

    return records


def dump_json(json_path, obj):
    with json_path.open('w') as fd:
        json.dump(obj, fd, indent=4, ensure_ascii=False)


def round_number(float_number):
    if float_number > 10:
        return Decimal(str(float_number)).quantize(Decimal('1'))
    elif float_number > 1:
        return Decimal(str(float_number)).quantize(Decimal('1.0'))
    elif float_number >= 0.1:
        return Decimal(str(float_number)).quantize(Decimal('1.00'))
    elif float_number == 0:
        return 0
    else:
        with localcontext() as ctx:
            ctx.prec = 2
            return Decimal(str(float_number)) + Decimal(0)
