SPECIAL_MUTATIONS = [
    '614G',
    '683G',
]


def include_similar_mutations(mutations):
    mut_list = []
    for mut in mutations:
        mut_list.append(mut)
        if not mut.startswith('S:'):
            continue
        for special in SPECIAL_MUTATIONS:
            if not mut.endswith(special):
                mut_list.append(mut + '+' + special)

    return include_mutations(mut_list)


def include_mutations(mutations):
    selector = ["s.iso_name = '{}'".format(m) for m in mutations]

    filter = """
        AND ({})
    """.format(' OR '.join(selector))

    return filter


def exclude_similar_mutations(mutations):
    mut_list = []
    for mut in mutations:
        mut_list.append(mut)
        if not mut.startswith('S:'):
            continue
        for special in SPECIAL_MUTATIONS:
            if not mut.endswith(special):
                mut_list.append(mut + '+' + special)

    return exclude_mutations(mut_list)


def exclude_mutations(mutations):
    selector = ["s.iso_name != '{}'".format(m) for m in mutations]

    filter = """
        AND ({})
    """.format(' AND '.join(selector))

    return filter
