from variant_filter import include_similar_mutations

from susceptibility import INDIVIDUAL_SAMPLE_SQL
from susceptibility import AGGREGATED_SUSC_VIEW_SQL
from variant.preset import CONTROL_VARIANTS_SQL
from variant.preset import get_iso_names_by_var_name


INDIVIDUAL_RESULTS_SQL = """
SELECT
    s.ref_name as ref_name,
    rx.cumulative_group as rx_name,
    s.control_iso_name as control,
    s.iso_name as iso_name,
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
WHERE
    s.control_iso_name IN {control_variants}
    {filters}
"""

INDIVIDUAL_CP_SQL = """
SELECT
    s.ref_name as ref_name,
    rx.cumulative_group as rx_name,
    rx.timing,
    rx.severity,
    s.control_iso_name as control,
    s.iso_name as iso_name,
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
WHERE
    s.control_iso_name IN {control_variants}
    {filters}
"""


INDIVIDUAL_VP_SQL = """
SELECT
    s.ref_name as ref_name,
    rx.vaccine_name as rx_name,
    rx.timing,
    rx.dosage,
    s.control_iso_name as control,
    s.iso_name as iso_name,
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
WHERE
    rx.vaccine_name IS NOT NULL AND
    s.control_iso_name IN {control_variants}
    {filters}
"""


AGGREGATED_RESULTS_SQL = """
SELECT
    s.ref_name as ref_name,
    rx.cumulative_group as rx_name,
    s.control_iso_name as control,
    s.iso_name as iso_name,
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
WHERE s.control_iso_name IN {control_variants}
    {filters}
"""


AGGREGATED_CP_SQL = """
SELECT
    s.ref_name as ref_name,
    rx.cumulative_group as rx_name,
    rx.timing,
    rx.severity,
    s.control_iso_name as control,
    s.iso_name as iso_name,
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
WHERE s.control_iso_name IN {control_variants}
    {filters}
"""

AGGREGATED_VP_SQL = """
SELECT
    s.ref_name as ref_name,
    rx.vaccine_name as rx_name,
    rx.timing,
    rx.dosage,
    s.control_iso_name as control,
    s.iso_name as iso_name,
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
WHERE
    rx.vaccine_name IS NOT NULL
    AND
    s.control_iso_name IN {control_variants}
    {filters}
"""


CP_FILTER = {
    'CP': {
        'rxtype': 'rx_conv_plasma',
        'cp_filters': [
            (
                "AND ("
                "      rx.infected_iso_name IN {}"
                "   OR rx.infected_iso_name IS NULL"
                "    )".format(CONTROL_VARIANTS_SQL)
            ),
        ]
    },
}

VP_FILTER = {
    'VP': {
        'rxtype': 'rx_vacc_plasma',
    },
}

RX_VP = """
(SELECT
    a.*,
    b.vaccine_type,
    b.priority
FROM
    rx_vacc_plasma AS a,
    vaccines AS b
ON
    a.vaccine_name == b.vaccine_name
)
"""


def get_iso_name_filter(var_name, selector='all'):

    iso_names = get_iso_names_by_var_name(var_name, selector)

    return "AND s.iso_name in ({})".format(iso_names)


VARIANTS = {
    'B.1.1.7': {
        'filter': [
            get_iso_name_filter('Alpha', 'spike')
        ]
    },
    'B.1.1.7 full genome': {
        'filter': [
            get_iso_name_filter('Alpha', 'genome')
        ]
    },
    'B.1.351': {
        'filter': [
            get_iso_name_filter('Beta', 'spike')
        ]
    },
    'B.1.351 full genome': {
        'filter': [
            get_iso_name_filter('Beta', 'genome')
        ]
    },
    'P.1': {
        'filter': [
            get_iso_name_filter('Gamma', 'spike')
        ]
    },
    'P.1 full genome': {
        'filter': [
            get_iso_name_filter('Gamma', 'genome')
        ]
    },
    'B.1.427/9': {
        'filter': [
            get_iso_name_filter('Epsilon')
        ]
    },
    'B.1.526': {
        'filter': [
            get_iso_name_filter('Iota')
        ]
    },
    'B.1.617.2': {
        'filter': [
            get_iso_name_filter('Delta')
        ]
    },
    'B.1.617.1': {
        'filter': [
            get_iso_name_filter('Kappa')
        ]
    },
}


MUTATIONS = {
    'N501Y': {
        'filter': [
            include_similar_mutations([
                'S:501Y',
            ])
        ]
    },
    '∆69/70': {
        'filter': [
            include_similar_mutations([
                'S:69del+70del',
            ])
        ]
    },
    '∆69/70 + N501Y': {
        'filter': [
            include_similar_mutations([
                'S:69del+70del+501Y',
            ])
        ]
    },
    '∆69/70 + N501Y + A570D': {
        'filter': [
            include_similar_mutations([
                'S:69del+70del+501Y+570D',
            ])
        ]
    },
    '∆69/70 + N501Y + Y453F': {
        'filter': [
            include_similar_mutations([
                'S:69del+70del+453F',
            ])
        ]
    },
    '∆144': {
        'filter': [
            include_similar_mutations([
                'S:144del',
            ])
        ]
    },
    'E484K': {
        'filter': [
            include_similar_mutations([
                'S:484K',
            ])
        ]
    },
    'E484K + N501Y': {
        'filter': [
            include_similar_mutations([
                'S:484K+501Y',
            ])
        ]
    },
    'Y453F': {
        'filter': [
            include_similar_mutations([
                'S:453F',
            ])
        ]
    },
    'L452R': {
        'filter': [
            include_similar_mutations([
                'S:452R',
            ])
        ]
    },
    'K417N': {
        'filter': [
            include_similar_mutations([
                'S:417N',
            ])
        ]
    },
    'K417N + E484K + N501Y': {
        'filter': [
            include_similar_mutations([
                'S:417N+484K+501Y',
                'B.1.351 RBD'
                ])
        ]
    },
    'N439K': {
        'filter': [
            include_similar_mutations([
                'S:439K',
            ])
        ]
    },
}


IGNORE_VACCINE_NAME = [
    'mRNA-1273.211',
    'mRNA-1273.351',
    'mRNA-1273 + mRNA-1273.351',
]

IGNORE_STUDY = [
    'Wall21b'
]
