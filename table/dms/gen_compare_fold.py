from collections import defaultdict
from operator import itemgetter
from preset import DATA_FILE_PATH
from preset import dump_csv


FOLD_SQL = """
SELECT
    susc.ref_name,
    susc.rx_name,
    rx.ab_name,
    susc.fold,
    susc.fold_cmp,
    susc.iso_name,
    mut.position,
    mut.amino_acid
FROM
    susc_results_50_wt_view susc,
    rx_mab_view rx,
    isolate_mutations_single_s_mut_view mut
WHERE
    susc.ref_name = rx.ref_name
    AND
    susc.rx_name = rx.rx_name
    AND
    susc.iso_name = mut.iso_name
    AND
    susc.potency_type = "IC50"
    AND
    rx.ab_name IN (SELECT ab_name FROM rx_dms_mab_view)
    AND
    susc.fold IS NOT NULL
GROUP BY
    susc.ref_name,
    susc.rx_name,
    susc.control_iso_name,
    susc.iso_name,
    susc.assay_name
;
"""


DMS_SQL = """
SELECT
    rx.ab_name,
    dms.position,
    dms.amino_acid,
    dms.escape_score
FROM
    dms_escape_results as dms,
    rx_dms_mab_view as rx
ON
    dms.ref_name = rx.ref_name AND
    dms.rx_name = rx.rx_name
WHERE
    position in ({positions})
"""


def gen_compare_fold(
        conn,
        save_path=DATA_FILE_PATH / 'dms' / 'fold_dms.csv'):

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
        positions=','.join(list([str(p) for p in positions]))
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
