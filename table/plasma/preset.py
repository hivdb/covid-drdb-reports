from variant_filter import include_mutations

from lookup_view import INDIVIDUAL_SAMPLE_SQL
from lookup_view import AGGREGATED_SUSC_VIEW_SQL


INDIVIDUAL_RESULTS_SQL = """
SELECT
    s.ref_name as ref_name,
    s.rx_name as rx_name,
    s.cumulative_count as sample_count,
    s.fold_cmp as fold_cmp,
    s.fold as fold
FROM
    (
""" + INDIVIDUAL_SAMPLE_SQL + """
) as s,
    {rxtype} AS rxtype
ON rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
WHERE s.control_variant_name IN ('Control', 'Wuhan', 'S:614G')
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
ON rxtype.ref_name = s.ref_name AND rxtype.rx_name = s.rx_name
WHERE s.control_variant_name IN ('Control', 'Wuhan', 'S:614G')
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
    'VP': {
        'rxtype': 'rx_immu_plasma',
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


EXCLUDE_STUDIES = {
    'Garcia-Beltran21': lambda x: x.startswith('BNT') or x.startswith('mRNA'),
    'Collier21': lambda x: x.startswith('CP'),
    'Widera21': lambda x: x.startswith('CP') or x.startswith('BNT'),
    'Brown21': lambda x: x.startswith('CP') or x.startswith('BNT'),
    'Supasa21': lambda x: x.startswith('CP')
                or x.startswith('BNT') or x.startswith('AZD'),
    'Edara21': lambda x: x == 'CP_acute',
}


EXCLUDE_PLASMA = [
    'Severe',
    'Mild',
]


RENAME_CP_EXECUTOR = {
    'Hoffmann21b': [
        (
            lambda x: x.startswith('ConvSerum') or x.startswith('ConvPlasma'),
            'CP',
        ),
    ],
    'Skelly21': [
        (
            lambda x: x.startswith('CP'),
            'CP'
        ),
        (
            lambda x: x.startswith('BNTsd'),
            'BNT162b2_one_dose'
        )
    ],
    'McCallum21b': [
        (
            lambda x: x.startswith('BNT_uninfected'),
            'BNT162b2'
        ),
        (
            lambda x: x.startswith('BNT_infected'),
            'BNT162b2_infected'
        ),
    ]
}


PLASMA_RENAME = {
    'BNT162b2_3W': 'BNT162b2_1M',
    'CP_5-33d': 'CP_1M',
    'CP_8M': 'CP_8M',
    'Mod_1M': 'mRNA-1273_1M',
    'Moderna_36d': 'mRNA-1273_1M',
    'Moderna_D43': 'mRNA-1273_1M',
    'NVV_1M': 'NVX-CoV2373_1M',
    'Pfizer_BNT162b2_D28': 'BNT162b2_1M',
    'Pfizer-BioNTech': 'BNT162b2',
    'Moderna': 'mRNA-1273',
    'BNT162b2_infected': 'BNT162b2',
    'CP_Mild': 'CP',
    "CP_Patient": 'CP',
    'CP_13': 'CP',
    'CP_29': 'CP',
    'CP_35': 'CP',
    'CP_37': 'CP',
    'CP_614G': 'CP',
    'CP_BNT162b2': 'BNT162b2',
    'CP_ModerateIgG': 'CP',
    'CP_StrongIgG': 'CP',
    'CP_weakIgG': 'CP',
    'Subject_B_d26': 'CP',
    'Subject_C_d32': 'CP',
    'Subject_G_d18': 'CP',
    'Subject_I_d26': 'CP',
    'CP_BNT162b2_1M': 'BNT162b2_1M',
    'BNT162b2_2W': 'BNT162b2_1M',
    'BNT162b2_4W': 'BNT162b2_1M',
    'CP_02-0014': 'CP_WT',
    'CP_02-0015': 'CP_WT',
    'CP_13-0013': 'CP_WT',
    'CP_13-0017': 'CP_WT',
    'CP_13-0033': 'CP_WT',
    'CP_13-0062': 'CP_WT',
    'BNT162b2_1M:01': 'BNT162b2',
    'BNT162b2_1M:02': 'BNT162b2',
    'BNT162b2_1M:03': 'BNT162b2',
    'BNT162b2_1M:04': 'BNT162b2',
    'BNT162b2_1M:05': 'BNT162b2',
    'BNT162b2_1M:06': 'BNT162b2',
    'BNT162b2_1M:07': 'BNT162b2',
    'BNT162b2_1M:08': 'BNT162b2',
    'BNT162b2_1M:09': 'BNT162b2',
    'BNT162b2_1M:10': 'BNT162b2',
    'BNT162b2_1M:11': 'BNT162b2',
    'BNT162b2_1M:12': 'BNT162b2',
    'BNT162b2_1M:13': 'BNT162b2',
    'BNT162b2_1M:14': 'BNT162b2',
    'BNT162b2_1M:15': 'BNT162b2',
    'CP_ID15': 'CP',
    'CP_ID18': 'CP',
    'CP_ID20': 'CP',
    'CP_ID22': 'CP',
    'CP_ID23': 'CP',
    'CP_ID24': 'CP',
    'CP_ID27': 'CP',
    'CP_ID33': 'CP',
    'CP_ID51': 'CP',
    'hCoV-2IG': 'IVIG',
}

PLASMA_POST_RENAME = {
    'BNT162b2_1M': 'BNT162b2',
    'CP_1M': 'CP',
    'CP_8M': 'CP(8M)',
    'CP_6M': 'CP(6M)',
    'NVX-CoV2373_1M': 'NVX-CoV',
    'mRNA-1273_1M': 'mRNA-1273',
    'CP_Titer_GT400': 'CP(High Titer)',
    'CP_Titer_LT400': 'CP(Low Titer)',
    'CP_WT': 'CP(WT)',
}
