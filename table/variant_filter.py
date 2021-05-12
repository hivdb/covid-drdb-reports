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


# def include_indiv_mutation(conn, mutation_str, include_614G=True):
#     cursor = conn.cursor()

#     SQL = """
#     SELECT iso_name, gene, position, amino_acid
#     FROM isolate_mutations
#     WHERE gene = "S"
#     AND iso_name IN (
#         SELECT iso_name
#         FROM isolate_mutations
#         WHERE gene = "S"
#         GROUP BY iso_name
#         HAVING count(1) in (1, 2))
#     """

#     final_iso_name = set()

#     for iso_name, gene, position, amino_acid in cursor.execute(SQL):
#         if position not in (501, 614):
#             if iso_name in final_iso_name:
#                 final_iso_name.remove(iso_name)
#         if position == 501 and amino_acid != 'Y':
#             if iso_name in final_iso_name:
#                 final_iso_name.remove(iso_name)
#         if position == 614 and amino_acid != 'G':
#             if iso_name in final_iso_name:
#                 final_iso_name.remove(iso_name)
