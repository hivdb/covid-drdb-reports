from variant.preset import ONE_MUT_VARIANT
from variant.preset import CONTROL_VARIANTS_SQL
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
    s.inhibition_pcnt != 90
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
    s.ordinal_number,
    s.assay
;
""".format(control_variants=CONTROL_VARIANTS_SQL)


DMS_SQL = """
SELECT
    d.rx_name,
    d.position,
    d.amino_acid,
    d.escape_score,
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
        save_path=DATA_FILE_PATH / 'summary_dms_missing_fold.csv'):

    cursor = conn.cursor()

    cursor.execute(SQL)

    mab2mutations = defaultdict(list)
    positions = set()

    for rec in cursor.fetchall():
        iso_name = rec['iso_name']
        if iso_name not in ONE_MUT_VARIANT.keys():
            continue

        mab = rec['rx_name']
        mab = MAB_RENAME.get(mab, mab)

        pos, aa = re.search(MUT_POS_AA, iso_name).groups()

        mab2mutations[mab].append((
            int(pos),
            aa,
        ))
        positions.add(pos)

    cursor.execute(DMS_SQL)

    results = []
    for rec in cursor.fetchall():
        mab = rec['rx_name']
        mab = MAB_RENAME.get(mab, mab)

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
