from collections import defaultdict
from json import dump
import json
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
    (
        (a.position = 339 AND a.amino_acid = 'D')
        OR
        (a.position = 371 AND a.amino_acid = 'L')
        OR
        (a.position = 371 AND a.amino_acid = 'F')
        OR
        (a.position = 373 AND a.amino_acid = 'P')
        OR
        (a.position = 375 AND a.amino_acid = 'F')
        OR
        (a.position = 376 AND a.amino_acid = 'A')
        OR
        (a.position = 405 AND a.amino_acid = 'N')
        OR
        (a.position = 408 AND a.amino_acid = 'S')
        OR
        (a.position = 417 AND a.amino_acid = 'N')
        OR
        (a.position = 440 AND a.amino_acid = 'K')
        OR
        (a.position = 446 AND a.amino_acid = 'S')
        OR
        (a.position = 477 AND a.amino_acid = 'N')
        OR
        (a.position = 478 AND a.amino_acid = 'K')
        OR
        (a.position = 484 AND a.amino_acid = 'A')
        OR
        (a.position = 493 AND a.amino_acid = 'R')
        OR
        (a.position = 496 AND a.amino_acid = 'S')
        OR
        (a.position = 498 AND a.amino_acid = 'R')
        OR
        (a.position = 501 AND a.amino_acid = 'Y')
        OR
        (a.position = 505 AND a.amino_acid = 'H')
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

        pos_group[position]['mut_aa'] = mut_aa
        pos_group[position]['ref_aa'] = ref_aa
        pos_group[position]['pos'] = position
        pos_group[position][ab_name] = escape_score

    records = []
    for _, ab_info in pos_group.items():
        records.append(ab_info)

    dump_csv(save_path, records)

    dump_json(json_save_path, records)
