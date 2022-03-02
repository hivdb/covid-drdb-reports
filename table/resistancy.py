import decimal
from decimal import Decimal
from decimal import localcontext

partial_cutoff = 5
resist_cutoff = 25


def get_susceptibility(fold_str):
    try:
        fold = float(fold_str)
        fold_cmp = '='
    except ValueError:
        fold_cmp = fold_str[0]
        fold = float(fold_str[1:])

    if fold < partial_cutoff:
        return 'susceptible'
    elif fold == partial_cutoff and fold_cmp in ['~', '<']:
        return 'susceptible'
    elif fold == partial_cutoff and fold_cmp in ['=', '>']:
        return 'partial-resistance'
    elif (fold > partial_cutoff and fold < resist_cutoff):
        return 'partial-resistance'
    elif fold == resist_cutoff and fold_cmp in ['~', '<']:
        return 'partial-resistance'
    elif fold == resist_cutoff and fold_cmp in ['=', '>']:
        return 'resistant'
    elif fold > resist_cutoff:
        return 'resistant'
    else:
        raise Exception(fold_str)


def is_susc(x):
    return x < partial_cutoff


def is_resistant(x):
    return x >= resist_cutoff


def is_partial_resistant(x):
    return x >= partial_cutoff and x < resist_cutoff


SUSCEPTIBLE_LEVEL_FILTER = """
    AND (
        (fold < 5
             OR (fold = 5 AND fold_cmp = '<')
             OR (fold = 5 AND fold_cmp = '~'))
        OR (resistance_level = 'susceptible')
        )
    AND (ineffective IS NULL)
"""

PARTIAL_RESISTANCE_LEVEL_FILTER = """
    AND (
        (
            (fold = 5 AND fold_cmp = '>')
            OR (fold = 5 AND fold_cmp = '=')
            OR (fold > 5 AND fold < 25)
            OR (fold = 25 AND fold_cmp = '~')
            OR (fold = 25 AND fold_cmp = '<'))
        OR (resistance_level = 'partial-resistance'))
    AND (ineffective IS NULL)
"""

RESISTANT_LEVLE_FILTER = """
    AND (
        (
            (fold = 25 AND fold_cmp = '>')
            OR (fold = 25 AND fold_cmp = '=')
            OR (fold > 25))
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
    elif number >= 10 and number < 1000:
        return Decimal(str(number)).quantize(Decimal('1'))
    # elif number == 1000:
    #     return 1000
    else:
        return '>1000'


def parse_fold(fold_str):
    try:
        fold = float(fold_str)
        # fold_cmp = '='
    except ValueError:
        # fold_cmp = fold_str[0]
        fold = float(fold_str[1:])
    return fold
