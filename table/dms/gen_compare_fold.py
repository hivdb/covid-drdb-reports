from variant.preset import ONE_MUT_VARIANT
from variant.preset import CONTROL_VARIANTS_SQL
from .preset import DMS_POSITIONS
from preset import DATA_FILE_PATH
from mab.preset import MAB_RENAME
from preset import dump_csv
import re
from collections import defaultdict
from operator import itemgetter


MUT_POS_AA = re.compile(r'(\d+)(\w)')


SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    s.fold,
    s.fold_cmp,
    s.iso_name
FROM
    susc_results as s,
    rx_antibodies as rx
ON
    s.ref_name = rx.ref_name
    AND s.rx_name = rx.rx_name
WHERE
    s.potency_type IN ('IC50', 'NT50')
    AND
    s.control_iso_name IN {control_variants}
    AND
    rx.ab_name IN (SELECT ab_name FROM rx_dms)
    AND
    s.fold IS NOT NULL
GROUP BY
    s.ref_name,
    s.rx_name,
    s.control_iso_name,
    s.iso_name,
    s.assay_name
;
""".format(control_variants=CONTROL_VARIANTS_SQL)


DMS_SQL = """
SELECT
    d.rx_name,
    position,
    amino_acid,
    escape_score
FROM
    dms_escape_results as d,
    rx_dms as rx
ON
    d.ref_name = rx.ref_name AND
    d.rx_name = rx.rx_name
WHERE
    position in ({positions})
"""


def gen_compare_fold(conn, save_path=DATA_FILE_PATH / 'summary_dms.csv'):

    cursor = conn.cursor()

    cursor.execute(SQL)

    indiv_mut_records = defaultdict(list)
    positions = set()

    for rec in cursor.fetchall():
        iso_name = rec['iso_name']
        if iso_name not in ONE_MUT_VARIANT.keys():
            continue

        ref_name = rec['ref_name']
        fold = rec['fold']
        fold_cmp = rec['fold_cmp']
        mab = rec['rx_name']
        mab = MAB_RENAME.get(mab, mab)

        pos, aa = re.search(MUT_POS_AA, iso_name).groups()

        indiv_mut_records[mab].append({
            'ref_name': ref_name,
            'mab': mab,
            'position': pos,
            'amino_acid': aa,
            'fold': fold,
            'fold_cmp': fold_cmp,
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
            ref_name = rec['ref_name']
            position = int(rec['position'])
            amino_acid = rec['amino_acid']
            fold = rec['fold']
            fold_cmp = rec['fold_cmp']

            pos_tree = eacape_tree.get(position)
            if not pos_tree:
                continue

            score = pos_tree.get(amino_acid)
            if not score:
                continue
            results.append({
                'ref_name': ref_name,
                'mab': mab,
                'position': position,
                'amino_acid': amino_acid,
                'fold': fold,
                'fold_cmp': fold_cmp,
                'score': score
            })

    results.sort(key=itemgetter('ref_name', 'mab'))

    dump_csv(save_path, results)
