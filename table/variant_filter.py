from variant.preset import ISONAME_MUTATIONS

SPECIAL_MUTATIONS = [
    '614G',
    '683G',
]


def include_similar_mutations(mutations):

    check_mut_list = set()
    for mut in mutations:
        gene, mutation_list = mut.split(':')
        check_mut_list.add(mutation_list)
        for m in SPECIAL_MUTATIONS:
            check_mut_list.add(mutation_list + '+' + m)

    iso_name_list = []
    for iso_name, iso_name_muts in ISONAME_MUTATIONS.items():
        # if iso_name_muts['non_s_mut_list']:
        #     continue
        if iso_name_muts['s_mut_str'] in check_mut_list:
            iso_name_list.append(iso_name)

    return include_mutations(iso_name_list)


def include_mutations(iso_name_list):
    iso_name_list = [
        "'{}'".format(iso_name)
        for iso_name in iso_name_list
    ]

    filter = """
        AND s.iso_name IN ({})
    """.format(','.join(iso_name_list))

    return filter


def exclude_similar_mutations(mutations):

    check_mut_list = set()
    for mut in mutations:
        gene, mutation_list = mut.split(':')
        check_mut_list.add(mutation_list)
        for m in SPECIAL_MUTATIONS:
            check_mut_list.add(mutation_list + '+' + m)

    iso_name_list = []
    for iso_name, iso_name_muts in ISONAME_MUTATIONS.items():
        # if iso_name_muts['non_s_mut_list']:
        #     continue
        if iso_name_muts['s_mut_str'] in check_mut_list:
            iso_name_list.append(iso_name)

    return exclude_mutations(iso_name_list)


def exclude_mutations(iso_name_list):
    iso_name_list = [
        "'{}'".format(iso_name)
        for iso_name in iso_name_list
    ]

    filter = """
        AND s.iso_name NOT IN ({})
    """.format(','.join(iso_name_list))

    return filter
