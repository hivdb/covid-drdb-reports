def include_mutations(mutations):
    selector = ["s.variant_name = '{}'".format(m) for m in mutations]

    filter = """
        AND ({})
    """.format(' OR '.join(selector))

    return filter


def exclude_mutations(mutations):
    selector = ["s.variant_name != '{}'".format(m) for m in mutations]

    filter = """
        AND ({})
    """.format(' AND '.join(selector))

    return filter


# def include_single_mutation(conn, mutation_str, include_614G=True):
#     cursor = conn.cursor()

#     SQL = """
#     SELECT variant_name, gene, position, amino_acid
#     FROM variant_mutations
#     WHERE gene = "S"
#     AND variant_name IN (
#         SELECT variant_name
#         FROM variant_mutations
#         WHERE gene = "S"
#         GROUP BY variant_name
#         HAVING count(1) in (1, 2))
#     """

#     final_variant_name = set()

#     for variant_name, gene, position, amino_acid in cursor.execute(SQL):
#         if position not in (501, 614):
#             if variant_name in final_variant_name:
#                 final_variant_name.remove(variant_name)
#         if position == 501 and amino_acid != 'Y':
#             if variant_name in final_variant_name:
#                 final_variant_name.remove(variant_name)
#         if position == 614 and amino_acid != 'G':
#             if variant_name in final_variant_name:
#                 final_variant_name.remove(variant_name)
