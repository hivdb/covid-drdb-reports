from collections import defaultdict
from operator import itemgetter
from preset import DATA_FILE_PATH
from preset import dump_csv


FOLD_SQL = """
SELECT
    susc.ref_name,
    susc.rx_name,
    susc.fold,
    susc.fold_cmp,
    susc.iso_name,
    rx.ab_name,
    mut.position,
    mut.amino_acid
FROM
    susc_results_50_wt_view as susc,
    rx_mab_view as rx,
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
    d.rx_name,
    d.position,
    d.amino_acid,
    d.escape_score,
    rx.ab_name,
    dms.ace2_binding,
    dms.expression,
    dms.ace2_contact
FROM
    dms_escape_results as d,
    rx_dms as rx,
    dms_ace2_binding as dms
ON
    d.ref_name = rx.ref_name AND
    d.rx_name = rx.rx_name AND
    d.position = dms.position AND
    d.amino_acid = dms.amino_acid
"""


def gen_missing_fold(
        conn,
        save_path=DATA_FILE_PATH / 'dms' / 'dms_missing_fold.csv'):

    cursor = conn.cursor()

    cursor.execute(FOLD_SQL)

    mab2mutations = defaultdict(list)
    positions = set()

    for rec in cursor.fetchall():
        pos = rec['position']
        aa = rec['amino_acid']

        mab = rec['ab_name']

        mab2mutations[mab].append((
            int(pos),
            aa,
        ))
        positions.add(pos)

    cursor.execute(DMS_SQL)

    results = []
    for rec in cursor.fetchall():
        mab = rec['ab_name']

        if mab not in mab2mutations.keys():
            continue

        mutations = mab2mutations[mab]

        position = int(rec['position'])
        amino_acid = rec['amino_acid']
        escape_score = rec['escape_score']
        if not escape_score:
            continue

        if float(escape_score) <= 0.1:
            continue

        if (position, amino_acid) in mutations:
            continue

        ace2_binding = rec['ace2_binding']
        expression = rec['expression']
        ace2_contact = rec['ace2_contact']

        results.append({
            'mab': mab,
            'position': position,
            'amino_acid': amino_acid,
            'score': escape_score,
            'ace2_binding': ace2_binding,
            'expression': expression,
            'ace2_contact': ace2_contact
        })

    results.sort(key=itemgetter('mab', 'position', 'amino_acid'))

    dump_csv(save_path, results)
