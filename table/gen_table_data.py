import sys
import sqlite3
from pathlib import Path

from variant.gen_table_key_variant import gen_table_key_variant
from variant.preset import get_grouped_variants
from variant.preset import get_spike_ref
from variant.gen_table_variant import gen_table_variant
from variant.gen_table_variant_mab import gen_table_variant_mab
from variant.gen_table_variant_vp import gen_table_variant_vp
from variant.gen_table_variant_cp import gen_table_variant_cp
from variant.gen_table_variant_aggre import gen_table_variant_aggre
from variant.gen_figure_variant_mab import gen_figure_variant_mab
from variant.gen_figure_variant_cp import gen_figure_variant_cp
from variant.gen_figure_variant_vp import gen_figure_variant_vp

from mab.gen_table_mab import gen_table_mab
from mab.gen_table_mab_variant import gen_table_mab_variant
from mab.gen_table_mab_muts import gen_table_mab_muts
from mab.gen_table_all_mab import gen_table_all_mab

from plasma.gen_table_cp_summary import gen_table_cp_summary
from plasma.gen_table_vp_summary import gen_table_vp_summary
from plasma.gen_table_cp_muts import gen_table_cp_muts
from plasma.gen_table_cp_variants import gen_table_cp_variants
from plasma.gen_table_vp_muts import gen_table_vp_muts
from plasma.gen_table_vp_variants import gen_table_vp_variants
from plasma.gen_vp_summary import gen_vp_summary
from plasma.gen_cp_summary import gen_cp_summary
from plasma.gen_plasma_figure import gen_plasma_figure

from study.gen_study import gen_study

from experiment.gen_null_fold import gen_null_fold
from experiment.gen_not_wildtype import gen_not_wildtype
from experiment.gen_exp import gen_exp

from dms.gen_compare_fold import gen_compare_fold
from dms.gen_issue_fold import gen_issue_fold


def gen_report(db_path):
    db_path = Path(db_path).resolve()
    print(db_path)
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row

    get_spike_ref(conn)

    get_grouped_variants(conn)

    gen_table_cp_muts(conn)
    gen_table_cp_variants(conn)
    gen_table_vp_muts(conn)
    gen_table_vp_variants(conn)
    gen_table_cp_summary()
    gen_table_vp_summary()
    gen_cp_summary(conn)
    gen_vp_summary(conn)
    gen_plasma_figure()

    gen_table_mab_variant(conn)
    gen_table_mab_muts(conn)
    gen_table_all_mab(conn)

    gen_table_key_variant(conn)
    gen_table_variant(conn)
    gen_table_variant_mab(conn)
    gen_table_variant_vp(conn)
    gen_table_variant_cp(conn)
    gen_table_variant_aggre(conn)
    gen_figure_variant_mab(conn)
    gen_figure_variant_vp(conn)
    gen_figure_variant_cp(conn)

    gen_study(conn)

    gen_null_fold(conn)
    gen_not_wildtype(conn)
    gen_exp(conn)

    gen_table_mab()

    gen_compare_fold(conn)
    gen_issue_fold()

    print('done')


if __name__ == '__main__':
    db_path = sys.argv[1]
    gen_report(db_path)
