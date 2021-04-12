def get_susceptibility(fold):
    try:
        fold = float(fold)
        fold_cmp = '='
    except ValueError:
        fold_cmp = fold[0]
        fold = float(fold[1:])

    if fold < 3:
        return 'susceptible'
    elif fold == 3 and fold_cmp in ['<', '=', '~']:
        return 'susceptible'
    elif fold == 3 and fold_cmp == '>':
        return 'partial-resistance'
    elif (fold > 3 and fold < 10):
        return 'partial-resistance'
    elif fold == 10 and fold_cmp in ['<', '=', '~']:
        return 'partial-resistance'
    else:
        return 'resistant'


def is_susc(x):
    return x <= 3


def is_resistant(x):
    return x >= 10


def is_partial_resistant(x):
    return x > 3 and x < 10


SUSCEPTIBLE_LEVEL_FILTER = """
    AND (
        (fold < 3
             OR (fold = 3 AND fold_cmp = '<')
             OR (fold = 3 AND fold_cmp = '=')
             OR (fold = 3 AND fold_cmp = '~'))
        OR (resistance_level = 'susceptible')
        )
    AND (ineffective IS NULL)
"""

PARTIAL_RESISTANCE_LEVEL_FILTER = """
    AND (
        (
            (fold = 3 AND fold_cmp = '>')
            OR (fold > 3 AND fold < 10)
            OR (fold = 10 AND fold_cmp = '=')
            OR (fold = 10 AND fold_cmp = '~')
            OR (fold = 10 AND fold_cmp = '<'))
        OR (resistance_level = 'partial-resistance'))
    AND (ineffective IS NULL)
"""

RESISTANT_LEVLE_FILTER = """
    AND (
        (
            (fold = 10 AND fold_cmp = '>')
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
