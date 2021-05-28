from variant_filter import include_mutations

from susceptibility import INDIVIDUAL_SAMPLE_SQL
from susceptibility import AGGREGATED_SUSC_VIEW_SQL
from variant.preset import CONTROL_VARIANTS_SQL
from .preset_cp import CP_VIEW

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


VARIANTS = {
    'B.1.1.7': {
        'filter': [
            "AND s.iso_name = 'B.1.1.7 Spike'",
        ]
    },
    'B.1.1.7 full genome': {
        'filter': [
            "AND s.iso_name = 'B.1.1.7 full genome'",
        ]
    },
    'B.1.351': {
        'filter': [
            "AND s.iso_name = 'B.1.351 Spike'",
        ]
    },
    'B.1.351 full genome': {
        'filter': [
            "AND s.iso_name = 'B.1.351 full genome'",
        ]
    },
    'P.1': {
        'filter': [
            "AND s.iso_name = 'P.1 Spike'",
        ]
    },
    'P.1 full genome': {
        'filter': [
            "AND s.iso_name = 'P.1 full genome'",
        ]
    },
    'B.1.427/9': {
        'filter': [
            "AND s.iso_name IN ("
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
    'B.1.617': {
        'filter': [
            include_mutations([
                'B.1.617 Spike',
                'B.1.617 full genome',
                'B.1.617.1_1675223',
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
    'BNT162b2': 'BNT162b2',
    'BNT_uninfected': 'BNT162b2',
    'BNT_W6': 'BNT162b2',
    'Mod': 'mRNA-1273',
    'MOD': 'mRNA-1273',
    'mRNA-1273': 'mRNA-1273',
    'mRNA_infected_1st': 'mRNA-1273',
    'mRNA_uninfected_1st': 'mRNA-1273',
    'mRNA_infected_2nd': 'mRNA-1273',
    'mRNA_uninfected_2nd': 'mRNA-1273',
    'Sputnik V': 'Sputnik V',
    'SPV': 'Sputnik V',
    'AZD': 'AZD1222',
    'AZD1222': 'AZD1222',
    'NVV': 'NVX-CoV',
    'NVX-CoV': 'NVX-CoV',
    'BBIBP-CorV': 'BBIBP-CorV',
    'CoronaVac': 'CoronaVac',
    'BBV152': 'BBV152',
    'MVC-COV1901': 'MVC-COV1901',
    'MVC': 'MVC-COV1901',
    'Mod_2nd': 'Mod',
    'BBIBP': 'BBIBP-CorV',
    'ZF2001': 'ZF2001',
}


VP_IGNORE = [
    ('Wu21b', 'Mod_211_mice'),
    ('Wu21b', 'Mod_351_mice'),
    ('Wu21b', 'Mod_mice_3rd'),
    ('Wu21b', 'Mod_mice'),
    ('Wu21b', 'Mod_mice_2nd'),
    ('Leier21*', 'BNT_infected'),
    ('McCallum21b', 'BNT_infected'),
    ('Collier21', 'BNT_1ds'),
    ('Supasa21*', 'AZD_14d'),
    ('Planas21', 'BNT_W2'),
    ('Planas21', 'BNT_W3'),
    ('Planas21', 'BNT_W4'),
    ('Planas21', 'BNT_W5'),
    ('Stankov21', 'BNT_1st_dose'),
    ('Pegu21', 'Mod_2nd_7M'),
    ('Pegu21*', 'Mod_2nd_7M'),
    ('Pegu21*', 'Mod_2nd_4M'),
    ('Pegu21*', 'Mod_1st'),
]
