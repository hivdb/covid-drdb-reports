from preset import DATA_FILE_PATH
from preset import dump_csv
from operator import itemgetter


SQL_TMPL = """
SELECT
    COUNT(DISTINCT rx.ref_name) num_ref,
    iso.domain
FROM
    susc_results as s,
    {rx_view} as rx,
    {iso_view} as iso
WHERE
    s.ref_name = rx.ref_name
    AND
    s.rx_name = rx.rx_name
    AND
    s.iso_name = iso.iso_name
GROUP BY
    iso.domain
"""


RX_VIEW = {
    'rx_conv_plasma': DATA_FILE_PATH / 'reference' / 'cp_domain.csv',
    'rx_vacc_plasma': DATA_FILE_PATH / 'reference' / 'vp_domain.csv',
    'rx_mab_view': DATA_FILE_PATH / 'reference' / 'mab_domain.csv',
}

ISO_TYPE_VIEW = {
    'single': 'isolate_mutations_single_s_mut_view',
    'combo': 'isolate_mutations_combo_s_mut_view',
}


def gen_ref_domain(conn):

    cursor = conn.cursor()

    for rx_view, save_path in RX_VIEW.items():

        results = []
        for iso_type, iso_view in ISO_TYPE_VIEW.items():

            sql = SQL_TMPL.format(
                rx_view=rx_view,
                iso_view=iso_view
                )

            cursor.execute(sql)
            for row in cursor.fetchall():

                domain = row['domain']
                num_ref = row['num_ref']
                results.append({
                    'variant_type': iso_type,
                    'domain': domain,
                    'num_ref': num_ref,
                })

        results.sort(key=itemgetter('variant_type', 'domain'))
        dump_csv(save_path, results)
