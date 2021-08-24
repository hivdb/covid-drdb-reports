import decimal
from decimal import Decimal
from decimal import localcontext


def get_susceptibility(fold_str):
    try:
        fold = float(fold_str)
        fold_cmp = '='
    except ValueError:
        fold_cmp = fold_str[0]
        fold = float(fold_str[1:])

    if fold < 3:
        return 'susceptible'
    elif fold == 3 and fold_cmp in ['~', '<']:
        return 'susceptible'
    elif fold == 3 and fold_cmp in ['=', '>']:
        return 'partial-resistance'
    elif (fold > 3 and fold < 10):
        return 'partial-resistance'
    elif fold == 10 and fold_cmp in ['~', '<']:
        return 'partial-resistance'
    elif fold == 10 and fold_cmp in ['=', '>']:
        return 'resistant'
    else:
        raise Exception(fold_str)


def is_susc(x):
    return x < 3


def is_resistant(x):
    return x >= 10


def is_partial_resistant(x):
    return x >= 3 and x < 10


SUSCEPTIBLE_LEVEL_FILTER = """
    AND (
        (fold < 3
             OR (fold = 3 AND fold_cmp = '<')
             OR (fold = 3 AND fold_cmp = '~'))
        OR (resistance_level = 'susceptible')
        )
    AND (ineffective IS NULL)
"""

PARTIAL_RESISTANCE_LEVEL_FILTER = """
    AND (
        (
            (fold = 3 AND fold_cmp = '>')
            OR (fold = 3 AND fold_cmp = '=')
            OR (fold > 3 AND fold < 10)
            OR (fold = 10 AND fold_cmp = '~')
            OR (fold = 10 AND fold_cmp = '<'))
        OR (resistance_level = 'partial-resistance'))
    AND (ineffective IS NULL)
"""

RESISTANT_LEVLE_FILTER = """
    AND (
        (
            (fold = 10 AND fold_cmp = '>')
            OR (fold = 10 AND fold_cmp = '=')
            OR (fold > 10))
        OR (resistance_level = 'resistant')
        OR (ineffective = 'experimental')
        )
"""

RESISTANCE_FILTER = {
    'susceptible': [SUSCEPTIBLE_LEVEL_FILTER],
    'partial': [PARTIAL_RESISTANCE_LEVEL_FILTER],
    'resistant': [RESISTANT_LEVLE_FILTER],
}


decimal.getcontext().rounding = decimal.ROUND_HALF_UP


def round_fold(number):
    if number < 10:
        number = float(number)
        if number.is_integer() and number != 0:
            return Decimal(str(number)).quantize(Decimal('1'))
        elif number >= 0.1:
            return Decimal(str(number)).quantize(Decimal('1.0'))
        else:
            return 0.1
            # with localcontext() as ctx:
            #     ctx.prec = 2
            #     return Decimal(str(number)) + Decimal(0)
    elif number >= 10 and number < 100:
        return Decimal(str(number)).quantize(Decimal('1'))
    # elif number == 100:
    #     return 100
    else:
        return '>100'
