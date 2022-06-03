from operator import itemgetter
from preset import DATA_FILE_PATH
from preset import dump_csv
from preset import dump_json
from sql import row2dict
from collections import defaultdict


def get_save_path(file_name, path, format):
    return path / '{}.{}'.format(file_name, format)


SQL = """
    SELECT
        DISTINCT
            iso.iso_name,
            iso.var_name,
            pat.pattern,
            iso.gisaid_id,
            iso.genbank_accn
    FROM
        isolate_mutations_variant_raw_s_mut_view pat,
        isolate_wildtype_view iso
    WHERE
        pat.iso_name = iso.iso_name
        AND
        pat.var_name = iso.var_name
        AND
        iso.iso_name IN
            (SELECT control_iso_name FROM susc_results_view)


    UNION ALL

    SELECT
        DISTINCT
            iso.iso_name,
            iso.var_name,
            '',
            iso.gisaid_id,
            iso.genbank_accn
    FROM
        isolate_wildtype_view iso
    WHERE
        NOT EXISTS (
            SELECT 1
            FROM
                isolate_mutations_variant_raw_s_mut_view pat
            WHERE
                pat.iso_name = iso.iso_name
                AND
                pat.var_name = iso.var_name
        )
        AND
        iso.iso_name IN
            (SELECT control_iso_name FROM susc_results_view)

;
"""


def gen_wildtype(
        conn,
        save_file_name='wildtype'
        ):

    cursor = conn.cursor()

    cursor.execute(SQL)

    records = row2dict(cursor.fetchall())

    csv_save_path = get_save_path(
        save_file_name, DATA_FILE_PATH / 'variant', 'csv')
    dump_csv(csv_save_path, records)

    var_pattern = defaultdict(set)
    for rec in records:
        var_name = rec['var_name']
        pattern = rec['pattern']
        var_pattern[var_name].add(pattern)

    records = []
    for var_name, pattern_list in var_pattern.items():
        for pattern in pattern_list:
            records.append({
                'var_name': var_name,
                'pattern': pattern
            })

    records.sort(key=itemgetter('var_name'))

    json_save_path = get_save_path(
        save_file_name, DATA_FILE_PATH / 'variant', 'json')
    dump_json(json_save_path, records)
