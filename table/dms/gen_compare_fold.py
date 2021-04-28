from variant.preset import INDIV_VARIANT
from variant.preset import CONTROL_VARIANTS_SQL
from .preset import DMS_POSITIONS
from preset import DATA_FILE_PATH
from mab.preset import MAB_RENAME
from preset import dump_csv
import re
from collections import defaultdict


MUT_POS_AA = re.compile('(\d+)(\w)')


SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.fold,
    s.fold_cmp,
    s.variant_name
FROM
    susc_results as s,
    rx_antibodies as rx
ON
    s.ref_name = rx.ref_name
    AND s.rx_name = rx.rx_name
WHERE
    s.control_variant_name IN {control_variants}
    AND
    rx.ab_name IN (SELECT ab_name FROM rx_dms)
    AND
    s.fold IS NOT NULL;
""".format(control_variants=CONTROL_VARIANTS_SQL)


DMS_SQL = """
SELECT
    rx_name,
    position,
    amino_acid,
    escape_score
FROM
    dms_escape_results
WHERE
    position in ({positions})
"""


def gen_compare_fold(conn, save_path=DATA_FILE_PATH / 'summary_dms.csv'):

    cursor = conn.cursor()

    cursor.execute(SQL)

    indiv_mut_records = defaultdict(list)
    positions = set()

    for rec in cursor.fetchall():
        variant_name = rec['variant_name']
        if variant_name not in INDIV_VARIANT.keys():
            continue

        fold = rec['fold']
        mab = rec['rx_name']
        mab = MAB_RENAME.get(mab, mab)

        pos, aa = re.search(MUT_POS_AA, variant_name).groups()

        indiv_mut_records[mab].append({
            'mab': mab,
            'position': pos,
            'amino_acid': aa,
            'fold': fold,
        })
        positions.add(pos)

    dms_sql = DMS_SQL.format(positions=','.join(list(positions)))

    cursor.execute(dms_sql)

    # ignored_rx = set()
    indiv_mut_escape = defaultdict(dict)
    for rec in cursor.fetchall():
        rx_name = rec['rx_name']
        rx_name = MAB_RENAME.get(rx_name, rx_name)
        if rx_name not in indiv_mut_records.keys():
            continue
        position = int(rec['position'])
        amino_acid = rec['amino_acid']
        escape_score = rec['escape_score']

        indiv_mut_escape[rx_name][position] = indiv_mut_escape[
            rx_name].get(position, {})

        indiv_mut_escape[rx_name][position][amino_acid] = escape_score

    results = []
    for mab, rx_list in indiv_mut_records.items():
        eacape_tree = indiv_mut_escape.get(mab)
        if not eacape_tree:
            continue

        for rec in rx_list:
            position = int(rec['position'])
            amino_acid = rec['amino_acid']
            fold = rec['fold']

            pos_tree = eacape_tree.get(position)
            if not pos_tree:
                continue

            score = pos_tree.get(amino_acid)
            if not score:
                continue
            results.append({
                'mab': mab,
                'position': position,
                'amino_acid': amino_acid,
                'fold': fold,
                'score': score
            })

    dump_csv(save_path, results)