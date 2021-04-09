import csv
from pathlib import Path
import decimal
from decimal import Decimal
import json

WS = Path(__file__).absolute().parent.parent
DATA_FILE_PATH = WS / 'susceptibility-data_files' / 'table'
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


decimal.getcontext().rounding = decimal.ROUND_HALF_UP


def round_number(number):
    if number < 10:
        number = float(number)
        if number.is_integer():
            return Decimal(str(number)).quantize(Decimal('1'))
        else:
            return Decimal(str(number)).quantize(Decimal('1.0'))
    elif number >= 10 and number < 100:
        return Decimal(str(number)).quantize(Decimal('1'))
    else:
        return '>100'
