from collections import defaultdict
from preset import dump_csv
from preset import dump_json
from preset import DATA_FILE_PATH

sql = """
SELECT
    a.gene,
    a.position position,
    a.amino_acid amino_acid,
    a.escape_score,
    b.ab_name,
    c.amino_acid ref_aa
FROM
    dms_escape_results a,
    dms_mab_view b,
    ref_amino_acid c
WHERE
    a.ref_name = b.ref_name
    AND
    a.rx_name = b.rx_name
    AND
    a.gene = 'S'
    AND
    a.position = c.position
    AND
    c.gene = 'S'
    AND
    b.availability IS NOT NULL
    AND
    EXISTS (
        SELECT 1
        FROM
            isolate_mutations d,
            isolates e
        WHERE
            d.iso_name = e.iso_name
            AND
            e.iso_name IN (
                "BA.1 Spike", "BA.1 Spike:493K",
                "BA.1 Spike:+346K", "BA.2 Spike")
            AND
            d.position = a.position
            AND
            d.amino_acid = a.amino_acid
            AND
            d.gene = 'S'
    )
;
"""


def gen_omicron(
    conn,
    save_path=DATA_FILE_PATH / 'variant' / 'omicron_dms.csv',
    json_save_path=DATA_FILE_PATH / 'variant' / 'omicron_dms.json'):

    cursor = conn.cursor()

    cursor.execute(sql)

    pos_group = defaultdict(dict)

    for row in cursor.fetchall():
        position = row['position']
        mut_aa = row['amino_acid']
        escape_score = row['escape_score']
        ab_name = row['ab_name']
        ref_aa = row['ref_aa']

        pos_group[(position, mut_aa)]['mut_aa'] = mut_aa
        pos_group[(position, mut_aa)]['ref_aa'] = ref_aa
        pos_group[(position, mut_aa)]['pos'] = position
        pos_group[(position, mut_aa)][ab_name] = escape_score

    records = []
    for _, ab_info in pos_group.items():
        records.append(ab_info)

    dump_csv(save_path, records)

    dump_json(json_save_path, records)
