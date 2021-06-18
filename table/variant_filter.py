def include_mutations(mutations):
    selector = ["s.iso_name = '{}'".format(m) for m in mutations]

    filter = """
        AND ({})
    """.format(' OR '.join(selector))

    return filter


def exclude_mutations(mutations):
    selector = ["s.iso_name != '{}'".format(m) for m in mutations]

    filter = """
        AND ({})
    """.format(' AND '.join(selector))

    return filter
