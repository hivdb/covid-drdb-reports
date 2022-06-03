import sys
import sqlite3
from pathlib import Path


def gen_report(db_path):
    db_path = Path(db_path).resolve()
    print(db_path)
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row

    from plasma.gen_figure_vp_variant import gen_figure_vp_variant
    from plasma.gen_cross_neutral import gen_cross_neutral

    from plasma.gen_table_cp_summary import gen_table_cp_summary
    from plasma.gen_table_vp_summary import gen_table_vp_summary
    from plasma.gen_table_cp_muts import gen_table_cp_muts
    from plasma.gen_table_cp_variants import gen_table_cp_variants
    from plasma.gen_table_vp_muts import gen_table_vp_muts
    from plasma.gen_table_vp_variants import gen_table_vp_variants

    from plasma.gen_vp_summary import gen_vp_summary
    from plasma.gen_vp_efficacy import gen_vp_efficacy
    from plasma.gen_vp_infection import gen_vp_infection
    from plasma.gen_vp_timing import gen_vp_timing
    from plasma.gen_vp_dosage import gen_vp_dosage
    from plasma.gen_vp_vaccine_species import gen_vp_vaccine_species
    from plasma.gen_vp_vaccine import gen_vp_vaccine

    from plasma.gen_cp_summary import gen_cp_summary
    from plasma.gen_cp_severity import gen_cp_severity
    from plasma.gen_cp_infection import gen_cp_infection
    from plasma.gen_cp_timing import gen_cp_timing

    from plasma.gen_figure_plasma_fold import gen_figure_plasma_fold
    from plasma.gen_figure_plasma_titer import gen_figure_plasma_titer
    from plasma.gen_figure_plasma_titer_fold import \
        gen_figure_plasma_titer_fold

    gen_table_cp_muts(conn)
    gen_table_cp_variants(conn)
    gen_table_vp_muts(conn)
    gen_table_vp_variants(conn)
    gen_table_cp_summary()
    gen_table_vp_summary()

    gen_vp_summary(conn)
    gen_vp_efficacy(conn)
    gen_vp_infection(conn)
    gen_vp_timing(conn)
    gen_vp_vaccine_species(conn)
    gen_vp_vaccine(conn)
    gen_vp_dosage(conn)

    gen_cp_summary(conn)
    gen_cp_severity(conn)
    gen_cp_infection(conn)
    gen_cp_timing(conn)

    gen_figure_plasma_fold(conn)
    gen_figure_plasma_titer(conn)
    gen_figure_plasma_titer_fold(conn)

    gen_figure_vp_variant(conn)
    gen_cross_neutral()

    from mab.gen_table_mab import gen_table_mab
    from mab.gen_table_mab_variant import gen_table_mab_variant
    from mab.gen_table_mab_muts import gen_table_mab_muts
    from mab.gen_table_all_mab import gen_table_all_mab

    from mab.gen_mab_synonyms import gen_mab_synonyms
    from mab.gen_mab_epitope_aligned import gen_mab_epitope_aligned

    gen_table_mab_variant(conn)
    gen_table_mab_muts(conn)
    gen_table_all_mab(conn)
    gen_table_mab()
    gen_mab_epitope_aligned(conn)
    gen_mab_synonyms(conn)

    from omicron.gen_omicron_ref_info import gen_omicron_ref_info
    from omicron.gen_omicron_titer_fold import gen_omicron_titer_fold
    from omicron.gen_omicron_assay import gen_omicron_assay
    from omicron.gen_omicron_single_mut import gen_omicron_single_mut
    from omicron.gen_table_mab_omicron_muts import gen_table_mab_omicron_muts
    from omicron.gen_table_mab_omicron import gen_table_mab_omicron
    from omicron.gen_ba_1_s_pattern import gen_ba_1_s_pattern
    from omicron.gen_omicron_wildtype_ic50 import gen_omicron_wildtype_ic50
    from omicron.gen_omicron_titer_fold_ba1_compare import \
        gen_omicron_titer_fold_ba1_compare
    from variant.gen_omicron import gen_omicron

    gen_omicron(conn)
    gen_omicron_titer_fold(conn)
    gen_omicron_ref_info(conn)
    gen_omicron_assay(conn)
    gen_ba_1_s_pattern(conn)
    gen_omicron_wildtype_ic50(conn)
    gen_omicron_single_mut(conn)
    gen_omicron_titer_fold_ba1_compare(conn)

    gen_table_mab_omicron_muts(conn)
    gen_table_mab_omicron()

    from variant.gen_table_variant_summary import gen_table_variant_summary
    from variant.gen_table_variant_mab import gen_table_variant_mab
    from variant.gen_table_variant_vp import gen_table_variant_vp
    from variant.gen_table_variant_cp import gen_table_variant_cp
    from variant.gen_table_variant_aggre import gen_table_variant_aggre
    from variant.gen_figure_variant_mab import gen_figure_variant_mab
    from variant.gen_figure_variant_cp import gen_figure_variant_cp
    from variant.gen_figure_variant_vp import gen_figure_variant_vp
    from variant.gen_table_single_mut_summary import \
        gen_table_single_mut_summary
    from variant.gen_wildtype import gen_wildtype

    gen_table_variant_summary(conn)
    gen_table_variant_mab(conn)
    gen_table_variant_vp(conn)
    gen_table_variant_cp(conn)
    gen_table_variant_aggre(conn)
    gen_figure_variant_mab(conn)
    gen_figure_variant_vp(conn)
    gen_figure_variant_cp(conn)
    gen_table_single_mut_summary(conn)
    gen_wildtype(conn)

    from reference.gen_ref_domain import gen_ref_domain
    from reference.gen_ref_info import gen_ref_info

    gen_ref_info(conn)
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
