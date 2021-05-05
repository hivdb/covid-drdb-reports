from collections import defaultdict

VARIANT_MUT_SQL = """
SELECT
    variant_name,
    position,
    amino_acid
FROM
    variant_mutations
WHERE gene = 'S';
"""

NO_MUT_SQL = """
SELECT
    variant_name
FROM
    virus_variants
WHERE variant_name
NOT IN (
    SELECT variant_name
    FROM variant_mutations
);
"""

CONTROL_VARIANTS_SQL = """
('Control', 'Wildtype', 'S:614G', 'WA1')
"""

IGNORE_MUTATION = [
    (614, 'G')
]

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
}

IGNORE_VARIANT_SYNYNOMS = [
    'Kemp21-d101',
    'WA-RML/d85',
    'WA-RML/d105',
    'B.1.427',
]

INDIV_VARIANT = {}
COMBO_VARIANT = {}
NO_MUT = []

NTD_DELETION = [
    (69, 70, '69-70∆'),
    (141, 145, '141-145∆'),
    (242, 244, '242-244∆'),
]


SPIKE_REF = {}
SPIKE_REF_SQL = """
SELECT * FROM 'ref_amino_acid' WHERE gene = 'S';
"""

DOMAINS = {
    'NTD': (1, 305),
    'RBD': (306, 534),
    'CTD': (535, 686),
    'S2': (687, 1273),
}


def get_spike_ref(conn):
    global SPIKE_REF

    cursor = conn.cursor()
    cursor.execute(SPIKE_REF_SQL)

    for rec in cursor.fetchall():
        position = rec['position']
        aa = rec['amino_acid']
        SPIKE_REF[position] = aa


def get_grouped_variants(conn):

    global INDIV_VARIANT
    global COMBO_VARIANT
    global NO_MUT

    cursor = conn.cursor()
    cursor.execute(VARIANT_MUT_SQL)

    variant_info = defaultdict(list)

    for rec in cursor.fetchall():
        variant = rec['variant_name']
        variant_info[variant].append(rec)

    uniq_variant_info = get_uniq_variant(variant_info)

    for mut_key, variant_info in uniq_variant_info.items():
        if variant_info['mut_count'] == 1:
            mutation = variant_info['mut_list'][0]
            INDIV_VARIANT[mutation['disp']] = mutation

            for name in variant_info['variant_names']:
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
            for name in variant_info['variant_names']:
                COMBO_VARIANT[name] = {
                    'disp': main_name,
                    'nickname': nickname,
                    }

    cursor.execute(NO_MUT_SQL)
    for rec in cursor.fetchall():
        NO_MUT.append(rec['variant_name'])

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

        uniq_variant['variant_names'] = uniq_variant.get('variant_names', [])
        uniq_variant['variant_names'].append(variant)

        uniq_variant['mut_list'] = mut_list

    for mut_key, uniq_variant in uniq_variant_info.items():
        mut_list = uniq_variant['mut_list']
        mut_list = merge_ntd_deletion(mut_list)
        uniq_variant['mut_list'] = mut_list
        uniq_variant['mut_count'] = len(mut_list)

    return uniq_variant_info


def get_main_name(variant):
    pass


def merge_ntd_deletion(mut_list):
    new_mut_list = []
    used_ntd_deletion = set()
    for mut in mut_list:
        if mut['aa'] != 'del':
            new_mut_list.append(mut)
            continue

        ref_aa = mut['ref_aa']
        position = mut['position']
        for ntd_start, ntd_stop, ntd_disp in NTD_DELETION:
            if ntd_disp in used_ntd_deletion:
                break
            if position >= ntd_start and position <= ntd_stop:
                new_mut_list.append({
                    'position': position,
                    'aa': 'del',
                    'ref_aa': ref_aa,
                    'disp': ntd_disp,
                    'domain': get_domain(position),
                })
                used_ntd_deletion.add(ntd_disp)

    return new_mut_list


def merge_spike_and_full_genome(variant_names):
    variant_names = [i.replace('Spike', '') for i in variant_names]
    variant_names = [i.replace('full genome', '') for i in variant_names]
    variant_names = [i.strip() for i in variant_names]

    variant_names = sorted(list(set(variant_names)))

    return variant_names


def get_combi_mutation_main_name(variant_info):
    mut_list = variant_info['mut_list']
    nickname = ''
    if len(mut_list) > 20:
        main_name = variant_info['variant_names'][0]
    else:
        main_name = ','.join([m['disp'] for m in mut_list])
        variant_names = variant_info['variant_names']
        for name in variant_names:
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
    for domain, range in DOMAINS.items():
        if position >= range[0] and position <= range[1]:
            return domain


def group_by_variant(records):
    indiv_records = defaultdict(list)
    combo_records = defaultdict(list)

    for rec in records:
        variant = rec['variant_name']
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
        variant_name = rec['variant_name']
        if variant_name in INDIV_VARIANT.keys():
            results.append(rec)
        elif variant_name == 'S:614G':
            results.append(rec)
        elif variant_name in NO_MUT:
            results.append(rec)
        elif variant_name in COMBO_VARIANT.keys():
            results.append(rec)

    return results
