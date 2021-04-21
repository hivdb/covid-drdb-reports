from variant_filter import include_mutations

from lookup_view import INDIVIDUAL_SAMPLE_SQL
from lookup_view import AGGREGATED_SUSC_VIEW_SQL

INDIVIDUAL_RESULTS_SQL = """
SELECT
    s.ref_name as ref_name,
    rx.cumulative_group as rx_name,
    s.control_variant_name as control,
    s.variant_name as variant_name,
    s.cumulative_count as sample_count,
    s.fold_cmp as fold_cmp,
    s.fold as fold
FROM
    (
""" + INDIVIDUAL_SAMPLE_SQL + """
) as s,
    {rxtype} AS rx
ON
    rx.ref_name = s.ref_name
    AND rx.rx_name = s.rx_name
WHERE s.control_variant_name IN {control_variants}
    {filters}
"""


AGGREGATED_RESULTS_SQL = """
SELECT
    s.ref_name as ref_name,
    rx.cumulative_group as rx_name,
    s.control_variant_name as control,
    s.variant_name as variant_name,
    s.cumulative_count as sample_count,
    s.fold_cmp as fold_cmp,
    s.fold as fold
FROM
    (
""" + AGGREGATED_SUSC_VIEW_SQL + """
) as s,
    {rxtype} AS rx
ON
    rx.ref_name = s.ref_name
    AND rx.rx_name = s.rx_name
WHERE s.control_variant_name IN {control_variants}
    {filters}
"""


SUBROWS = {
    'CP': {
        'rxtype': 'rx_conv_plasma',
        'cp_filters': [
            (
                "AND ("
                "      rx.infection IN ('Wildtype', 'S:614G')"
                "   OR rx.infection IS NULL"
                "    )"
            ),
        ]
    },
    'VP': {
        'rxtype': 'rx_vacc_plasma',
    },
}


VARIANTS = {
    'B.1.1.7': {
        'filter': [
            "AND s.variant_name = 'B.1.1.7 Spike'",
        ]
    },
    'B.1.1.7 full genome': {
        'filter': [
            "AND s.variant_name = 'B.1.1.7 full genome'",
        ]
    },
    'B.1.351': {
        'filter': [
            "AND s.variant_name = 'B.1.351 Spike'",
        ]
    },
    'B.1.351 full genome': {
        'filter': [
            "AND s.variant_name = 'B.1.351 full genome'",
        ]
    },
    'P.1': {
        'filter': [
            "AND s.variant_name = 'P.1 Spike'",
        ]
    },
    'P.1 full genome': {
        'filter': [
            "AND s.variant_name = 'P.1 full genome'",
        ]
    },
    'B.1.427/9': {
        'filter': [
            "AND s.variant_name IN ("
            "    'B.1.427 full genome',"
            "    'B.1.429 full genome',"
            "    'B.1.429 Spike')",
        ]
    },
    'B.1.526': {
        'filter': [
            include_mutations([
                'B.1.526 Spike',
                'B.1.526 full genome',
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


VP_RENAME = {
    'BNT': 'BNT162b2',
    'MOD': 'mRNA-1273',
    'SPV': 'Sputnik V',
    'AZD': 'AZD1222',
    'NVV': 'NVX-CoV',
}
