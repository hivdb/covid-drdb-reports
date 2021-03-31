from variant_filter import include_mutations

from lookup_view import INDIVIDUAL_SAMPLE_SQL
from lookup_view import AGGREGATED_SUSC_VIEW_SQL


INDIVIDUAL_RESULTS_SQL = """
SELECT
    s.ref_name as ref_name,
    s.rx_name as rx_name,
    s.cumulative_count as sample_count,
    s.fold as fold
FROM
    (
""" + INDIVIDUAL_SAMPLE_SQL + """
) as s,
    {rxtype} AS rxtype
    WHERE rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
    AND s.control_variant_name IN ('Control', 'Wuhan', 'S:614G')
    -- AND s.ineffective IS NULL
    {filters}
"""


AGGREGATED_RESULTS_SQL = """
SELECT
    s.ref_name as ref_name,
    s.rx_name as rx_name,
    s.cumulative_count as sample_count,
    s.fold_cmp as fold_cmp,
    s.fold as fold
FROM
    (
""" + AGGREGATED_SUSC_VIEW_SQL + """
) as s,
    {rxtype} AS rxtype
    WHERE rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
    AND s.control_variant_name IN ('Control', 'Wuhan', 'S:614G')
    -- AND s.ineffective IS NULL
    {filters}
    GROUP BY s.ref_name, s.rx_name;
"""


SUBROWS = {
    'CP': {
        'rxtype': 'rx_conv_plasma',
        'cp_filters': [
            (
                "AND ("
                "      rxtype.infection IN ('S:614G')"
                "   OR rxtype.infection IS NULL"
                "    )"
            ),
        ]
    },
    'IP': {
        'rxtype': 'rx_immu_plasma',
    },
}


VARIANTS = {
    'B.1.1.7': {
        'filter': [
            "AND s.variant_name = 'B.1.1.7 Spike'",
        ]
    },
    'B.1.1.7 authentic': {
        'filter': [
            "AND s.variant_name = 'B.1.1.7 authentic'",
        ]
    },
    'B.1.351': {
        'filter': [
            "AND s.variant_name = 'B.1.351 Spike'",
        ]
    },
    'B.1.351 authentic': {
        'filter': [
            "AND s.variant_name = 'B.1.351 authentic'",
        ]
    },
    'P.1': {
        'filter': [
            "AND s.variant_name = 'P.1 Spike'",
        ]
    },
    'P.1 authentic': {
        'filter': [
            "AND s.variant_name = 'P.1 authentic'",
        ]
    },
    'B.1.427/9': {
        'filter': [
            "AND s.variant_name IN ("
            "    'B.1.427 authentic',"
            "    'B.1.429 authentic',"
            "    'B.1.429 Spike')",
        ]
    },
    'B.1.526': {
        'filter': [
            include_mutations([
                'B.1.526 Spike',
                'B.1.526 authentic',
            ])
        ]
    },
}


MUTATIONS = {
    'N501Y': {
        'filter': [
            include_mutations([
                'S:501Y',
                'S:501Y+614G'])
        ]
    },
    '∆69/70': {
        'filter': [
            include_mutations([
                'S:69del+70del',
                'S:69del+70del+614G'])
        ]
    },
    '∆69/70 + N501Y': {
        'filter': [
            include_mutations([
                'S:69del+70del+501Y',
                'S:69del+70del+501Y+614G'])
        ]
    },
    '∆69/70 + N501Y + A570D': {
        'filter': [
            include_mutations([
                'S:69del+70del+501Y+570D',
                'S:69del+70del+501Y+570D+614G'])
        ]
    },
    '∆69/70 + N501Y + Y453F': {
        'filter': [
            include_mutations([
                'S:69del+70del+453F',
                'S:69del+70del+453F+614G'])
        ]
    },
    '∆144': {
        'filter': [
            include_mutations([
                'S:144del',
                'S:144del+614G'])
        ]
    },
    'E484K': {
        'filter': [
            include_mutations([
                'S:484K',
                'S:484K+614G'])
        ]
    },
    'E484K + N501Y': {
        'filter': [
            include_mutations([
                'S:484K+501Y',
                'S:484K+501Y+614G'])
        ]
    },
    'Y453F': {
        'filter': [
            include_mutations([
                'S:453F',
                'S:453F+614G'])
        ]
    },
    'L452R': {
        'filter': [
            include_mutations([
                'S:452R',
                'S:452R+614G'])
        ]
    },
    'K417N': {
        'filter': [
            include_mutations([
                'S:417N',
                'S:417N+614G'])
        ]
    },
    'K417N + E484K + N501Y': {
        'filter': [
            include_mutations([
                'S:417N+484K+501Y',
                'S:417N+484K+501Y+614G',
                'B.1.351 RBD'
                ])
        ]
    },
    'N439K': {
        'filter': [
            include_mutations([
                'S:439K',
                'S:439K+614G'])
        ]
    },
}
