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

IGNORE_MUTATION = [
    (614, 'G')
]

IGNORE_VARIANT_SYNYNOMS = [
    'Kemp21-d101',
    'WA-RML/d85',
    'WA-RML/d105',
    'B.1.427',
]

INDIV_VARIANT = {}
MULTI_VARIANT = {}
NO_MUT = []


def get_grouped_variants(conn):

    global INDIV_VARIANT
    global MULTI_VARIANT
    global NO_MUT

    cursor = conn.cursor()
    cursor.execute(VARIANT_MUT_SQL)

    variant_info = defaultdict(list)

    for rec in cursor.fetchall():
        variant = rec['variant_name']
        variant_info[variant].append(rec)

    uniq_variant = get_uniq_variant(variant_info)

    for i in uniq_variant:
        main_name = i['main_name']
        if i['mut_count'] == 1:
            for name in i['variant_names']:
                INDIV_VARIANT[name] = main_name
        else:
            for name in i['variant_names']:
                MULTI_VARIANT[name] = main_name

    cursor.execute(NO_MUT_SQL)
    for rec in cursor.fetchall():
        NO_MUT.append(rec['variant_name'])


def get_uniq_variant(variant_info):

    uniq_variant_info = defaultdict(dict)

    for variant, info_list in variant_info.items():
        mut_list = []
        for rec in info_list:
            position = rec['position']
            aa = rec['amino_acid']
            if (position, aa) in IGNORE_MUTATION:
                continue
            mut_list.append((position, aa))

        mut_list = sorted(mut_list, key=lambda x: x[0])
        mut_count = len(mut_list)
        mut_key = ','.join(
            ['{}{}'.format(mut[0], mut[1]) for mut in mut_list])

        uniq_variant = uniq_variant_info[mut_key]

        uniq_variant['mut_count'] = mut_count

        uniq_variant['variant_names'] = uniq_variant.get('variant_names', [])
        uniq_variant['variant_names'].append(variant)

        # uniq_variant['mut_list'] = uniq_variant.get('mut_list', [])
        # uniq_variant['mut_list'].append(mut_list)

    for mut_key, uniq_variant in uniq_variant_info.items():
        variant_names = uniq_variant['variant_names']
        variant_names = [name.replace('+614G', '') for name in variant_names]
        variant_names = [name.replace('614G+', '') for name in variant_names]
        variant_names = sorted(variant_names)
        variant_names = merge_spike_and_full_genome(variant_names)

        if len(variant_names) == 1:
            uniq_variant['main_name'] = variant_names[0]
        else:
            variant_names = [
                i for i in variant_names if i not in IGNORE_VARIANT_SYNYNOMS]
            if len(variant_names) == 1:
                uniq_variant['main_name'] = variant_names[0]
            else:
                variant_names = [
                    i for i in variant_names if not i.startswith('S:')]
                uniq_variant['main_name'] = variant_names[0]

    return list(uniq_variant_info.values())


def merge_spike_and_full_genome(variant_names):
    variant_names = [i.replace('Spike', '') for i in variant_names]
    variant_names = [i.replace('full genome', '') for i in variant_names]
    variant_names = [i.strip() for i in variant_names]

    variant_names = sorted(list(set(variant_names)))

    return variant_names
