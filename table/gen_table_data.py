import sys
import sqlite3
from pathlib import Path


def gen_report(db_path):
    db_path = Path(db_path).resolve()
    print(db_path)
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row

    from variant.preset import get_spike_ref
    get_spike_ref(conn)

    from variant.preset import get_grouped_variants
    get_grouped_variants(conn)

    from variant.preset import gen_iso_name2var_name
    gen_iso_name2var_name(conn)

    from variant.preset import load_isoname_mutations
    load_isoname_mutations(conn)

    # from plasma.preset import CP_VIEW
    # print(CP_VIEW)
    from plasma.gen_table_cp_summary import gen_table_cp_summary
    from plasma.gen_table_vp_summary import gen_table_vp_summary
    from plasma.gen_table_cp_muts import gen_table_cp_muts
    from plasma.gen_table_cp_variants import gen_table_cp_variants
    from plasma.gen_table_vp_muts import gen_table_vp_muts
    from plasma.gen_table_vp_variants import gen_table_vp_variants
    from plasma.gen_vp_summary import gen_vp_summary
    from plasma.gen_cp_summary import gen_cp_summary
    from plasma.gen_figure_plasma_fold import gen_figure_plasma_fold
    from plasma.gen_figure_plasma_titer import gen_figure_plasma_titer
    from plasma.gen_figure_plasma_titer_fold import \
        gen_figure_plasma_titer_fold
    from plasma.gen_vp_efficacy import gen_vp_efficacy

    gen_table_cp_muts(conn)
    gen_table_cp_variants(conn)
    gen_table_vp_muts(conn)
    gen_table_vp_variants(conn)
    gen_table_cp_summary()
    gen_table_vp_summary()
    gen_cp_summary(conn)
    gen_vp_summary(conn)
    gen_vp_efficacy(conn)
    gen_figure_plasma_fold(conn)
    gen_figure_plasma_titer(conn)
    gen_figure_plasma_titer_fold(conn)

    from mab.gen_table_mab import gen_table_mab
    from mab.gen_table_mab_variant import gen_table_mab_variant
    from mab.gen_table_mab_muts import gen_table_mab_muts
    from mab.gen_table_all_mab import gen_table_all_mab

    gen_table_mab_variant(conn)
    gen_table_mab_muts(conn)
    gen_table_all_mab(conn)
    gen_table_mab()

    from variant.gen_table_variant_summary import gen_table_variant_summary
    from variant.gen_table_variant_mab import gen_table_variant_mab
    from variant.gen_table_variant_vp import gen_table_variant_vp
    from variant.gen_table_variant_cp import gen_table_variant_cp
    from variant.gen_table_variant_aggre import gen_table_variant_aggre
    from variant.gen_figure_variant_mab import gen_figure_variant_mab
    from variant.gen_figure_variant_cp import gen_figure_variant_cp
    from variant.gen_figure_variant_vp import gen_figure_variant_vp

    gen_table_variant_summary(conn)
    gen_table_variant_mab(conn)
    gen_table_variant_vp(conn)
    gen_table_variant_cp(conn)
    gen_table_variant_aggre(conn)
    gen_figure_variant_mab(conn)
    gen_figure_variant_vp(conn)
    gen_figure_variant_cp(conn)

    from reference.gen_ref_domain import gen_ref_domain

    gen_ref_domain(conn)

    from fold.gen_null_fold import gen_null_fold
    from fold.gen_non_wt_control import gen_non_wt_control
    from fold.gen_assay import gen_assay
    from fold.gen_control import gen_control

    gen_null_fold(conn)
    gen_non_wt_control(conn)
    gen_assay(conn)
    gen_control(conn)

    from dms.gen_compare_fold import gen_compare_fold
    from dms.gen_missing_fold import gen_missing_fold
    from dms.gen_issue_fold import gen_issue_fold

    gen_compare_fold(conn)
    gen_missing_fold(conn)
    gen_issue_fold()

    print('done')


if __name__ == '__main__':
    db_path = sys.argv[1]
    gen_report(db_path)
