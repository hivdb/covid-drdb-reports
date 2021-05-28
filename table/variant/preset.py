from collections import defaultdict

DOMAINS = {
    'NTD': (1, 305),
    'RBD': (306, 534),
    'CTD': (535, 686),
    'S2': (687, 1273),
}

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


CONTROL_VARIANTS_SQL = None

D614G_ONLY_MUTATION = """
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

def gen_control_variants(conn):
    global CONTROL_VARIANTS_SQL

    cursor = conn.cursor()
    cursor.execute(D614G_ONLY_MUTATION)
    isolates = [r for r in cursor.fetchall()]

    cursor.execute(NO_S_MUTATION)
    isolates += [r for r in cursor.fetchall()]

    # print('Control isolates', [
    #     (r['iso_name'], r['var_name']) for r in isolates])

    iso_names = [r['iso_name'] for r in isolates]

    CONTROL_VARIANTS_SQL = '({})'.format(
        ', '.join(["'{}'".format(i) for i in iso_names])
    )


INDIV_VARIANT = {}
COMBO_VARIANT = {}
NO_MUT = []

VARIANT_MUT_SQL = """
SELECT
    iso_name,
    position,
    amino_acid
FROM
    isolate_mutations
WHERE gene = 'S';
"""

NO_MUT_SQL = """
SELECT
    iso_name
FROM
    isolates
WHERE iso_name
NOT IN (
    SELECT iso_name
    FROM isolate_mutations
);
"""

IGNORE_VARIANTS = [
    'SARS-CoV Spike',
    'L455Y,F456L',
    'T470N,E471V,I472P',
    'E324G,S325D',
    'K444T,V445S,G446T',
    'WIV1 Spike',
]

VARIANT_NICKNAMES = {
    '69-70∆,N501Y': 'B.1.1.7 (variation)',
    '69-70∆,N501Y,P681H': 'B.1.1.7 (variation)',
    '69-70∆,N501Y,A570D': 'B.1.1.7 (variation)',
    'K417T,E484K,N501Y': 'P.1 (variation)',
    'K417T,E484K': 'P.1 (variation)',
    'K417N,N501Y': 'B.1.351 (variation)',
    'K417N,E484K': 'B.1.351 (variation)',
    'B.1.427': 'B.1.427/9',
    'B.1.1.7 S1': 'B.1.1.7 (variation)',
    'B.1.1.7 + 484K (variation)': 'B.1.1.7 + E484K',
    'B.1.1.7 :-144del+145del (variation)': 'B.1.1.7',
    'B.1.351 :-18F-215G (variation)': 'B.1.351 (variation)',
    'B.1.351 :-18F-242del-243del-244del-246I (variation)':
        'B.1.351 (variation)',
    'B.1.351 :-18F-242del-243del-244del-246I-417N (variation)':
        'B.1.351 (variation)',
    'B.1.351 :-18F-244del-246I (variation)': 'B.1.351 (variation)',
    'B.1.351 :-242del-243del-244del-246I (variation)': 'B.1.351 (variation)',
    'B.1.351 :-246I (variation)': 'B.1.351 (variation)',
    'B.1.351 RBD': 'B.1.351 (variation)',
    'B.1.526 :-484K+477N (variation)': 'B.1.526 (variation)',
    'B.1.526 v.2': 'B.1.526 (variation)',
}

NTD_DELETION_GROUP_RULE = {
    (69, 70): '69-70∆',
    (69, 69): '69-70∆',
    (139, 145): '141-145∆',
    (140, 140): '141-145∆',
    (141, 145): '141-145∆',
    (141, 144): '141-145∆',
    (144, 144): '141-145∆',
    (145, 145): '141-145∆',
    (145, 146): '145-146∆',
    (141, 143): '141-145∆',
    (147, 147): '147∆',
    (241, 243): '242-244∆',
    (242, 244): '242-244∆',
    (244, 244): '242-244∆',
    (244, 247): '244-247∆',
    (245, 246): '244-247∆',
    (254, 257): '254-257∆',
    (72, 78): '72-78∆',
    (150, 152): '150-152∆',
}

IGNORE_MUTATION = [
    (614, 'G')
]

IGNORE_VARIANT_SYNYNOMS = [
    'Kemp21-d101',
    'WA-RML/d85',
    'WA-RML/d105',
    'B.1.427',
]


def get_grouped_variants(conn):

    global INDIV_VARIANT
    global COMBO_VARIANT
    global NO_MUT

    cursor = conn.cursor()
    cursor.execute(VARIANT_MUT_SQL)

    variant_info = defaultdict(list)

    for rec in cursor.fetchall():
        variant = rec['iso_name']
        variant_info[variant].append(rec)

    uniq_variant_info = get_uniq_variant(variant_info)

    for mut_key, variant_info in uniq_variant_info.items():
        if variant_info['mut_count'] == 1:
            mutation = variant_info['mut_list'][0]
            INDIV_VARIANT[mutation['disp']] = mutation

            for name in variant_info['iso_names']:
                INDIV_VARIANT[name] = mutation
        else:
            main_name, nickname = get_combi_mutation_main_name(variant_info)
            if main_name in IGNORE_VARIANTS:
                continue

            nickname = VARIANT_NICKNAMES.get(
                main_name,
                VARIANT_NICKNAMES.get(nickname, nickname))
            COMBO_VARIANT[main_name] = {
                'disp': main_name,
                'nickname': nickname,
                }
            for name in variant_info['iso_names']:
                COMBO_VARIANT[name] = {
                    'disp': main_name,
                    'nickname': nickname,
                    }

    cursor.execute(NO_MUT_SQL)
    for rec in cursor.fetchall():
        NO_MUT.append(rec['iso_name'])

    # from pprint import pprint
    # pprint(list(COMBO_VARIANT.keys()))


def get_uniq_variant(variant_info):

    uniq_variant_info = defaultdict(dict)

    for variant, info_list in variant_info.items():
        mut_list = []
        for rec in info_list:
            position = rec['position']
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
            continue

        mut_key = ','.join(
            [mut['disp'] for mut in mut_list])

        uniq_variant = uniq_variant_info[mut_key]
        uniq_variant['mut_count'] = mut_count

        uniq_variant['iso_names'] = uniq_variant.get('iso_names', [])
        uniq_variant['iso_names'].append(variant)

        uniq_variant['mut_list'] = mut_list

    for mut_key, uniq_variant in uniq_variant_info.items():
        mut_list = uniq_variant['mut_list']
        mut_list = merge_ntd_deletion(mut_list)
        uniq_variant['mut_list'] = mut_list
        uniq_variant['mut_count'] = len(mut_list)

    return uniq_variant_info


def get_combi_mutation_main_name(variant_info):
    mut_list = variant_info['mut_list']
    nickname = ''
    if len(mut_list) > 20:
        main_name = variant_info['iso_names'][0]
    else:
        main_name = ','.join([m['disp'] for m in mut_list])
        iso_names = variant_info['iso_names']
        for name in iso_names:
            if nickname:
                break
            if name.startswith('S:'):
                continue

            nickname = name.replace('Spike', '').strip()
            nickname = nickname.replace('full genome', '').strip()
            if 'S:' in nickname:
                nickname = nickname.replace('S:', ' ')
                nickname += ' (variation)'
            if ':-' in nickname:
                nickname += ' (variation)'

    return main_name, nickname


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
    iso_names = [i.replace('full genome', '') for i in iso_names]
    iso_names = [i.strip() for i in iso_names]

    iso_names = sorted(list(set(iso_names)))

    return iso_names


def group_by_variant(records):
    indiv_records = defaultdict(list)
    combo_records = defaultdict(list)

    for rec in records:
        variant = rec['iso_name']
        main_name = INDIV_VARIANT.get(variant)
        if main_name:
            disp_name = main_name['disp']
            indiv_records[disp_name].append(rec)
        else:
            main_name = COMBO_VARIANT.get(variant)
            if not main_name:
                continue
            disp_name = main_name['disp']
            combo_records[disp_name].append(rec)

    return indiv_records, combo_records


def filter_by_variant(records):
    results = []
    for rec in records:
        iso_name = rec['iso_name']
        if iso_name in INDIV_VARIANT.keys():
            results.append(rec)
        elif iso_name == 'S:614G':
            results.append(rec)
        elif iso_name in NO_MUT:
            results.append(rec)
        elif iso_name in COMBO_VARIANT.keys():
            results.append(rec)

    return results
