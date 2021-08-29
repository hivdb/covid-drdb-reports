from variant.preset import CONTROL_VARIANTS_SQL
from variant.preset import SINGLE_S_MUTATION_ISOLATES
from preset import DATA_FILE_PATH
from mab.preset import RX_MAB
from mab.preset import RX_MAB_DMS
from preset import dump_csv
import re
from collections import defaultdict
from operator import itemgetter


MUT_POS_AA = re.compile(r'(\d+)(\w)')


FOLD_SQL = """
SELECT
    s.ref_name,
    s.rx_name,
    rx.ab_name,
    s.fold,
    s.fold_cmp,
    s.iso_name,
    iso.position,
    iso.amino_acid
FROM
    susc_results as s,
    ({rx_antibodies}) as rx,
    ({isolate_mutations}) as iso
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.iso_name = iso.iso_name
    AND
    s.potency_type LIKE "IC%"
    AND
    s.control_iso_name IN {control_variants}
    AND
    rx.ab_name IN (SELECT ab_name FROM ({rx_mab_dms}))
    AND
    s.fold IS NOT NULL
GROUP BY
    s.ref_name,
    s.rx_name,
    s.control_iso_name,
    s.iso_name,
    s.assay_name
;
""".format(
    control_variants=CONTROL_VARIANTS_SQL,
    rx_antibodies=RX_MAB,
    rx_mab_dms=RX_MAB_DMS,
    isolate_mutations=SINGLE_S_MUTATION_ISOLATES
    )


DMS_SQL = """
SELECT
    rx.ab_name,
    dms.position,
    dms.amino_acid,
    dms.escape_score
FROM
    dms_escape_results as dms,
    ({rx_mab_dms}) as rx
ON
    dms.ref_name = rx.ref_name AND
    dms.rx_name = rx.rx_name
WHERE
    position in ({positions})
"""


def gen_compare_fold(conn, save_path=DATA_FILE_PATH / 'summary_dms.csv'):

    cursor = conn.cursor()

    cursor.execute(FOLD_SQL)

    indiv_mut_records = defaultdict(list)
    positions = set()

    for rec in cursor.fetchall():
        ref_name = rec['ref_name']
        fold = rec['fold']
        fold_cmp = rec['fold_cmp']
        mab = rec['ab_name']
        pos = int(rec['position'])
        aa = rec['amino_acid']

        indiv_mut_records[mab].append({
            'ref_name': ref_name,
            'mab': mab,
            'position': pos,
            'amino_acid': aa,
            'fold': fold,
            'fold_cmp': fold_cmp,
        })
        positions.add(pos)

    dms_sql = DMS_SQL.format(
        positions=','.join(list(
            [str(p) for p in positions])),
        rx_mab_dms=RX_MAB_DMS
    )

    cursor.execute(dms_sql)

    # ignored_rx = set()
    indiv_mut_escape = defaultdict(dict)
    for rec in cursor.fetchall():
        mab = rec['ab_name']
        if mab not in indiv_mut_records.keys():
            continue

        position = int(rec['position'])
        amino_acid = rec['amino_acid']
        escape_score = rec['escape_score']

        if mab not in indiv_mut_escape:
            indiv_mut_escape[mab] = defaultdict(dict)

        if position not in indiv_mut_escape[mab]:
            indiv_mut_escape[mab][position] = {}

        indiv_mut_escape[mab][position][amino_acid] = escape_score

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
