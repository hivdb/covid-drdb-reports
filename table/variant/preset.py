from collections import defaultdict

DOMAINS = {
    'NTD': (1, 305),
    'RBD': (306, 534),
    'CTD': (535, 686),
    'S2': (687, 1273),
}

ISOLATE_MUTATIONS_WITH_REF = """
SELECT
    a.iso_name,
    a.gene,
    b.amino_acid as ref,
    a.position,
    a.amino_acid
FROM
    isolate_mutations a,
    ref_amino_acid b
WHERE
    a.gene = b.gene
    AND
    a.position = b.position
"""


CONTROL_VARIANTS_SQL = """
SELECT
    a.iso_name
FROM
    isolates a,
    variants b
WHERE
    a.var_name = b.var_name
    AND
    b.as_wildtype
"""

D614G_ONLY_ISOLATE = """
SELECT
    b.iso_name,
    b.var_name
FROM
    (
    SELECT
        iso_name,
        GROUP_CONCAT(mutation, ',') as mutations
    FROM (
        SELECT
            iso_name,
            position || amino_acid as mutation
        FROM
            isolate_mutations
        WHERE
            gene = 'S'
        )
    GROUP BY
        iso_name
    )
    AS a,
    isolates
    AS b
ON
    a.iso_name = b.iso_name
WHERE
    mutations = '614G'
;
"""

SPIKE_ONLY_ISOLATES = """
SELECT
    *
FROM
    ({isolate_mutations}) a
WHERE NOT EXISTS (
    SELECT
        1
    FROM
        ({isolate_mutations}) b
    WHERE
        gene != 'S'
        AND
        b.iso_name = a.iso_name
    )
""".format(
    isolate_mutations=ISOLATE_MUTATIONS_WITH_REF)

SPIKE_ONLY_MUT_ISOLATES = """
SELECT
    *
FROM
    ({spike_only_isolates})
WHERE
    (position, amino_acid) != (614, 'G')
    AND
    (position, amino_acid) != (683, 'G')
""".format(
    spike_only_isolates=SPIKE_ONLY_ISOLATES
)


SINGLE_S_MUTATION_ISOLATES = """
SELECT
    ref || position || amino_acid single_mut_name,
    iso_name,
    ref,
    position,
    amino_acid
FROM
    ({spike_only_mut_isolates}) a
WHERE
    EXISTS (
        SELECT
            *
        FROM
            ({spike_only_mut_isolates}) b
        WHERE
            a.iso_name = b.iso_name
        GROUP BY
            b.iso_name
        HAVING
            COUNT(1) = 1
    )
""".format(
    spike_only_mut_isolates=SPIKE_ONLY_MUT_ISOLATES)


NO_S_MUTATION = """
SELECT
    iso_name,
    var_name
FROM
    isolates
WHERE
    iso_name != 'Unknown'
    AND
    iso_name NOT IN
    (
        SELECT
            distinct iso_name
        FROM
            isolate_mutations
        WHERE
            gene = 'S'
    )
"""

WT_MUTATIONS = """
SELECT
    *
FROM
    isolates
WHERE
    var_name in ('A', 'B', 'B.1')
    OR iso_name in ('S:614G', 'S:683G')
;
"""


ISO_NAME_SQL = """
SELECT
    iso_name,
    var_name
FROM
    isolates
"""

ISO_NAME2VAR_NAME = {}
VAR_NAME2ISO_NAMES = defaultdict(list)


def gen_iso_name2var_name(conn):
    global ISO_NAME2VAR_NAME
    global VAR_NAME2ISO_NAMES

    cursor = conn.cursor()
    cursor.execute(ISO_NAME_SQL)

    for rec in cursor.fetchall():
        iso_name = rec['iso_name']
        var_name = rec['var_name']

        if not var_name:
            continue
        ISO_NAME2VAR_NAME[iso_name] = var_name

        VAR_NAME2ISO_NAMES[var_name].append(iso_name)


def get_iso_names_by_var_name(var_name, selector='all'):
    iso_names = VAR_NAME2ISO_NAMES.get(var_name, [])

    if not iso_names:
        return ''

    spike_only_iso_names = [
        name for name, info in ISONAME_MUTATIONS.items()
        if not len(info['non_s_mut_list'])
        ]

    if selector == 'all':
        pass
    elif selector == 'spike':
        iso_names = [s for s in iso_names if s in spike_only_iso_names]
    elif selector == 'genome':
        iso_names = [s for s in iso_names if s not in spike_only_iso_names]

    return ','.join(["'{}'".format(s) for s in iso_names])


ONE_MUT_VARIANT = {}
COMBO_MUT_VARIANT = {}
NO_S_MUT_VARIANT = []

VARIANT_MUT_SQL = """
SELECT
    iso_muts.iso_name,
    iso.var_name,
    iso_muts.position,
    iso_muts.amino_acid
FROM
    isolate_mutations AS iso_muts,
    isolates AS iso
ON
    iso_muts.iso_name = iso.iso_name
WHERE
    gene = 'S';
"""

IGNORE_VARIANTS = [
    'SARS-CoV Spike',
    'L455Y,F456L',
    'T470N,E471V,I472P',
    'E324G,S325D',
    'K444T,V445S,G446T',
    'WIV1 Spike',
]

IGNORE_VARIANT_SYNYNOMS = [
    'Kemp21-d101',
    'WA-RML/d85',
    'WA-RML/d105',
    'B.1.427',
]

NTD_DELETION_GROUP_RULE = {
    (69, 70): '69-70∆',
    (69, 69): '69-70∆',
    (63, 75): '63-75∆',
    (139, 145): '141-145∆',
    (140, 140): '141-145∆',
    (141, 141): '141-145∆',
    (141, 145): '141-145∆',
    (141, 144): '141-145∆',
    (144, 144): '141-145∆',
    (144, 145): '141-145∆',
    (145, 145): '141-145∆',
    (145, 146): '145-146∆',
    (141, 143): '141-145∆',
    (147, 147): '147∆',
    (210, 210): '210∆',
    (241, 241): '242-244∆',
    (241, 243): '242-244∆',
    (242, 243): '242-244∆',
    (242, 244): '242-244∆',
    (244, 244): '242-244∆',
    (244, 247): '244-247∆',
    (245, 246): '244-247∆',
    (246, 248): '244-247∆',
    (254, 257): '254-257∆',
    (72, 78): '72-78∆',
    (150, 152): '150-152∆',
    (156, 157): '157-158∆',
    (157, 158): '157-158∆',
    (156, 158): '157-158∆',
    (246, 252): '247-253∆',
    (247, 253): '247-253∆',
}

IGNORE_MUTATION = [
    (614, 'G')
]


def get_grouped_variants(conn):

    global ONE_MUT_VARIANT
    global COMBO_MUT_VARIANT
    global NO_S_MUT_VARIANT

    cursor = conn.cursor()
    cursor.execute(VARIANT_MUT_SQL)

    isolate_group = defaultdict(list)

    for rec in cursor.fetchall():
        iso_name = rec['iso_name']
        isolate_group[iso_name].append(rec)

    uniq_isolate_list = get_uniq_isolate_pattern_list(isolate_group)

    for mut_key, isolate_info in uniq_isolate_list.items():
        if isolate_info['mut_count'] == 1:
            mutation = isolate_info['mut_list'][0]
            ONE_MUT_VARIANT[mutation['disp']] = mutation

            for name in isolate_info['iso_names']:
                ONE_MUT_VARIANT[name] = mutation
        else:
            if len(isolate_info['var_names']) > 1:
                print(mut_key)
                print('iso_name', isolate_info['iso_names'])
                print('var_name', isolate_info['var_names'])

            main_name = get_combined_mutation_main_name(isolate_info)
            if main_name in IGNORE_VARIANTS:
                continue

            var_name = ','.join(isolate_info['var_names'])

            COMBO_MUT_VARIANT[main_name] = {
                'disp': main_name,
                'varname': var_name,
                }
            for name in isolate_info['iso_names']:
                COMBO_MUT_VARIANT[name] = {
                    'disp': main_name,
                    'varname': var_name,
                    }

    cursor.execute(NO_S_MUTATION)
    for rec in cursor.fetchall():
        NO_S_MUT_VARIANT.append(rec['iso_name'])

    # from pprint import pprint
    # pprint(list(COMBO_MUT_VARIANT.keys()))


SPIKE_REF = {}
SPIKE_REF_SQL = """
SELECT
    *
FROM
    ref_amino_acid
WHERE
    gene = 'S'
;
"""


def get_spike_ref(conn):
    global SPIKE_REF

    cursor = conn.cursor()
    cursor.execute(SPIKE_REF_SQL)

    for rec in cursor.fetchall():
        position = rec['position']
        aa = rec['amino_acid']
        SPIKE_REF[position] = aa


def get_uniq_isolate_pattern_list(isolate_info):

    uniq_isolate_list = defaultdict(dict)

    for iso_name, records in isolate_info.items():
        mut_list = []

        uniq_var_name_list = set()

        for rec in records:
            var_name = rec['var_name'] or ''
            uniq_var_name_list.add(var_name)

            position = int(rec['position'])
            aa = rec['amino_acid']
            ref_aa = SPIKE_REF[position]
            if (position, aa) in IGNORE_MUTATION:
                continue

            mut_list.append({
                'ref_aa': ref_aa,
                'position': position,
                'aa': aa,
                'disp': '{}{}{}'.format(ref_aa, position, aa),
                'domain': get_domain(position)
            })

        mut_list = sorted(mut_list, key=lambda x: x['position'])

        mut_count = len(mut_list)
        if mut_count == 0:
            # print('No mutations', iso_name)
            continue

        uniq_mutation_list = ','.join(
            [mut['disp'] for mut in mut_list])

        uniq_isolate = uniq_isolate_list[uniq_mutation_list]
        uniq_isolate['mut_count'] = mut_count

        uniq_isolate['iso_names'] = uniq_isolate.get('iso_names', [])
        uniq_isolate['iso_names'].append(iso_name)

        uniq_isolate['var_names'] = uniq_isolate.get('var_names', set())
        uniq_isolate['var_names'] |= uniq_var_name_list

        uniq_isolate['mut_list'] = mut_list

    for uniq_mutation_list, uniq_isolate in uniq_isolate_list.items():
        mut_list = uniq_isolate['mut_list']

        mut_list = merge_ntd_deletion(mut_list)

        uniq_isolate['mut_list'] = mut_list
        uniq_isolate['mut_count'] = len(mut_list)

    return uniq_isolate_list


def get_combined_mutation_main_name(variant_info):
    mut_list = variant_info['mut_list']
    main_name = ','.join([m['disp'] for m in mut_list])
    return main_name


def get_domain(position):
    position = int(position)
    for domain, range in DOMAINS.items():
        if position >= range[0] and position <= range[1]:
            return domain


def merge_ntd_deletion(mut_list):
    no_del_mut_list = [
        m for m in mut_list
        if (m['aa'] != 'del' and m['domain'] == 'NTD')
        or (m['domain'] != 'NTD')
        ]
    del_mut_list = [
        m for m in mut_list
        if (m['aa'] == 'del' and m['domain'] == 'NTD')
        ]

    del_mut_list.sort(key=lambda x: int(x['position']))

    consequitive_del_group = []
    prev_position = -1
    current_del_group = []
    for mut in del_mut_list:
        position = int(mut['position'])
        if position != prev_position + 1:
            if current_del_group:
                consequitive_del_group.append(current_del_group)

            current_del_group = []

        prev_position = position
        current_del_group.append(mut)

    if current_del_group:
        consequitive_del_group.append(current_del_group)

    grouped_del_mut_list = []

    for del_group in consequitive_del_group:
        start_pos = del_group[0]['position']
        stop_pos = del_group[-1]['position']

        rule_key = (start_pos, stop_pos)

        ntd_disp = NTD_DELETION_GROUP_RULE[rule_key]
        ref_aas = ''.join([m['ref_aa'] for m in del_group])
        grouped_del_mut_list.append({
            'position': start_pos,
            'aa': 'del',
            'ref_aa': ref_aas,
            'disp': ntd_disp,
            'domain': get_domain(start_pos),
        })

    new_mut_list = grouped_del_mut_list + no_del_mut_list

    new_mut_list.sort(key=lambda x: int(x['position']))

    return new_mut_list


def merge_spike_and_full_genome(iso_names):
    iso_names = [i.replace('Spike', '') for i in iso_names]
    iso_names = [i[:i.find('full genome')] for i in iso_names]
    iso_names = [i.strip() for i in iso_names]

    iso_names = sorted(list(set(iso_names)))

    return iso_names


def group_by_variant(records):
    indiv_records = defaultdict(list)
    combo_records = defaultdict(list)

    for rec in records:
        variant = rec['iso_name']
        main_name = ONE_MUT_VARIANT.get(variant)
        if main_name:
            disp_name = main_name['disp']
            indiv_records[disp_name].append(rec)
        else:
            main_name = COMBO_MUT_VARIANT.get(variant)
            if not main_name:
                continue
            disp_name = main_name['disp']
            combo_records[disp_name].append(rec)

    return indiv_records, combo_records


def filter_by_variant(records):
    results = []
    for rec in records:
        iso_name = rec['iso_name']
        if iso_name in ONE_MUT_VARIANT.keys():
            results.append(rec)
        elif iso_name == 'S:614G':
            results.append(rec)
        elif iso_name in NO_S_MUT_VARIANT:
            results.append(rec)
        elif iso_name in COMBO_MUT_VARIANT.keys():
            results.append(rec)

    return results


ISONAME_MUTATIONS = {}

ISO_MUT_QUERY = """
SELECT iso_name, gene, position, amino_acid
FROM isolate_mutations;
"""


def load_isoname_mutations(conn):
    global ISONAME_MUTATIONS

    cursor = conn.cursor()
    cursor.execute(ISO_MUT_QUERY)

    iso_name_mut = defaultdict(list)

    for rec in cursor.fetchall():
        iso_name = rec['iso_name']
        iso_name_mut[iso_name].append({
            'gene': rec['gene'],
            'pos_aa': rec['position'],
            'mut_aa': rec['amino_acid'],
        })

    for iso_name, mut_list in iso_name_mut.items():
        s_mut_list = [m for m in mut_list if m['gene'] == 'S']
        non_s_mut_list = [m for m in mut_list if m['gene'] != 'S']
        s_mut_str = '+'.join(
            [
                '{}{}'.format(
                    m['pos_aa'],
                    m['mut_aa']
                )
                for m in s_mut_list
            ]
        )
        ISONAME_MUTATIONS[iso_name] = {}
        ISONAME_MUTATIONS[iso_name]['s_mut_str'] = s_mut_str
        ISONAME_MUTATIONS[iso_name]['non_s_mut_list'] = non_s_mut_list
