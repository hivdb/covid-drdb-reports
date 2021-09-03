import sqlite3
from pathlib import Path
from preset import load_csv, dump_csv
from preset import load_json
from mab import RX_MAB_DMS
from mab import RX_MAB_DRDB
from mab import RX_SINGLE_MAB_DMS_SQL
from variants import SINGLE_SPIKE_MUTATION
from variants import WT_VARIANTS
from vp import RX_FULLY_VACCINE
from statistics import median
from statistics import quantiles
from preset import round_number


POSITION = 449
AA = ''


AA_NAMES = 'ACDEFGHIKLMNPQRSTVWY'

DB_PATH = Path(__file__).parent.absolute() / 'covid-drdb-latest.db'
OUTBREAK_PATH = Path(__file__).parent.absolute() / 'variants-mutations.csv'
OUTBREAK_COUNT_PATH = Path(__file__).parent.absolute() / 'variants-count.json'


def load_variant_count():
    result = {}
    for i in load_json(OUTBREAK_COUNT_PATH)['results']:
        name = i['name'].upper()
        total = i['total_count']
        result[name] = total
    return result


def run_query(sql_stat):
    conn = sqlite3.connect(str(DB_PATH))
    conn.row_factory = sqlite3.Row

    cursor = conn.cursor()
    cursor.execute(sql_stat)

    result = []
    for row in cursor.fetchall():
        rec = {}
        for key in row.keys():
            rec[key] = row[key]
        result.append(rec)
    return result


def query_dms_by_mutation(pos, aa):
    print('\n==== DMS ====')
    DMS_BINDING_SQL = """
    SELECT
        ace2_binding,
        expression
    FROM
        dms_ace2_binding
    WHERE
        position = {pos}
        AND
        amino_acid = '{aa}'
    """.format(pos=pos, aa=aa)

    for i in run_query(DMS_BINDING_SQL):
        for k, v in i.items():
            print(k, ':', v)

    DMS_SCORE_SQL = """
    SELECT
        a.ab_name,
        b.escape_score
    FROM
        ({single_dms_mab}) a,
        dms_escape_results b
    WHERE
        a.rx_name = b.rx_name
        AND
        a.availability IS NOT NULL
        AND
        position = {pos}
        AND
        amino_acid = '{aa}'
    """.format(
        single_dms_mab=RX_SINGLE_MAB_DMS_SQL,
        pos=pos,
        aa=aa)

    for i in run_query(DMS_SCORE_SQL):
        show_values = []
        for v in i.values():
            show_values.append(str(v))
        print(':\t'.join(show_values))


def query_invitro_by_mutation(pos, aa):
    print('\n==== In Vitro Selection ====')
    INVITRO_SQL = """
    SELECT
        b.ab_name,
        group_concat(a.ref_name, ',') ref_names
    FROM
        invitro_selection_results a,
        ({rx_mab_drdb}) b
    WHERE
        a.ref_name = b.ref_name
        AND
        a.rx_name = b.rx_name
        AND
        b.availability IS NOT NULL
        AND
        a.position = {pos}
        AND
        a.amino_acid = '{aa}'
    GROUP BY
        b.ab_name
    """.format(
        rx_mab_drdb=RX_MAB_DRDB,
        pos=pos,
        aa=aa
    )
    for i in run_query(INVITRO_SQL):
        show_values = []
        for v in i.values():
            show_values.append(str(v))
        print(':\t'.join(show_values))


def _format_print_rx_info(records):
    rx_info = {}
    for i in records:
        rx_name = i['main_name']
        fold = i['fold']
        num_results = i['num_results']
        if rx_name not in rx_info:
            rx_info[rx_name] = {
                'fold_list': [],
                'num_results': 0
            }

        rx_info[rx_name]['fold_list'].append(fold)
        rx_info[rx_name]['num_results'] += num_results

    print('  |  '.join(['Rx', 'median', 'min=>max', 'iqr', 'num_results']))

    for rx_name, _info in rx_info.items():
        med = median(_info['fold_list'])
        med = round_number(med)
        max_value = max(_info['fold_list'])
        min_value = min(_info['fold_list'])
        range = '{}=>{}'.format(
            round_number(min_value),
            round_number(max_value)
            )
        if len(_info['fold_list']) >= 4:
            iqr = quantiles(_info['fold_list'])
            iqr = '[{}, {}]'.format(
                round_number(iqr[0]),
                round_number(iqr[-1])
                )
        else:
            iqr = '-'
        print('  |  '.join([
            str(i) for i in
            [rx_name, med, range, iqr, _info['num_results']]]))


def query_invivo_by_mutation(pos, aa):
    print('\n==== mAb invitro (Clinical trials) ====')
    MAB_EUA_SQL = """
    SELECT
        rx.ab_name main_name,
        s.fold_cmp,
        s.fold,
        s.iso_name,
        sum(s.cumulative_count) num_results
    FROM
        susc_results s,
        ({rx_mab}) rx
    WHERE
        s.ref_name = rx.ref_name
        AND
        s.rx_name = rx.rx_name
        AND
        s.control_iso_name IN ({control_iso_sql})
        AND
        rx.availability IS NOT NULL
        AND
        s.iso_name IN (
            SELECT
                iso_name
            FROM
                ({single_mutation_isolates})
            WHERE
                position = {pos}
                AND
                amino_acid = '{aa}'
        )
    GROUP BY
        s.ref_name,
        s.rx_name
    """.format(
        rx_mab=RX_MAB_DRDB,
        single_mutation_isolates=SINGLE_SPIKE_MUTATION,
        control_iso_sql=WT_VARIANTS,
        pos=pos,
        aa=aa
    )

    _format_print_rx_info(run_query(MAB_EUA_SQL))

    print('\n==== mAb invitro (other) ====')
    MAB_NON_EUA_SQL = """
    SELECT
        rx.ab_name main_name,
        s.fold_cmp,
        s.fold,
        s.iso_name,
        sum(s.cumulative_count) num_results
    FROM
        susc_results s,
        ({rx_mab}) rx
    WHERE
        s.ref_name = rx.ref_name
        AND
        s.rx_name = rx.rx_name
        AND
        s.control_iso_name IN ({control_iso_sql})
        AND
        rx.availability IS NULL
        AND
        s.iso_name IN (
            SELECT
                iso_name
            FROM
                ({single_mutation_isolates})
            WHERE
                position = {pos}
                AND
                amino_acid = '{aa}'
        )
    GROUP BY
        s.ref_name,
        s.rx_name
    """.format(
        rx_mab=RX_MAB_DRDB,
        single_mutation_isolates=SINGLE_SPIKE_MUTATION,
        control_iso_sql=WT_VARIANTS,
        pos=pos,
        aa=aa
    )

    _format_print_rx_info(run_query(MAB_NON_EUA_SQL))

    print('\n==== CP invitro ====')
    CP_SQL = """
    SELECT
        'CP' main_name,
        s.fold_cmp,
        s.fold,
        s.iso_name,
        sum(s.cumulative_count) num_results
    FROM
        susc_results s,
        rx_conv_plasma rx
    WHERE
        s.ref_name = rx.ref_name
        AND
        s.rx_name = rx.rx_name
        AND
        s.control_iso_name IN ({control_iso_sql})
        AND
        rx.infected_iso_name IN ({control_iso_sql})
        AND
        s.iso_name IN (
            SELECT
                iso_name
            FROM
                ({single_mutation_isolates})
            WHERE
                position = {pos}
                AND
                amino_acid = '{aa}'
        )
    GROUP BY
        s.ref_name,
        s.rx_name
    """.format(
        single_mutation_isolates=SINGLE_SPIKE_MUTATION,
        control_iso_sql=WT_VARIANTS,
        pos=pos,
        aa=aa
    )

    _format_print_rx_info(run_query(CP_SQL))

    print('\n==== VP invitro ====')
    VP_SQL = """
    SELECT
        'VP' main_name,
        s.fold_cmp,
        s.fold,
        s.iso_name,
        sum(s.cumulative_count) num_results
    FROM
        susc_results s,
        ({rx_vp}) rx
    WHERE
        s.ref_name = rx.ref_name
        AND
        s.rx_name = rx.rx_name
        AND
        s.control_iso_name IN ({control_iso_sql})
        AND
        s.iso_name IN (
            SELECT
                iso_name
            FROM
                ({single_mutation_isolates})
            WHERE
                position = {pos}
                AND
                amino_acid = '{aa}'
        )
    GROUP BY
        s.ref_name,
        s.rx_name
    """.format(
        rx_vp=RX_FULLY_VACCINE,
        single_mutation_isolates=SINGLE_SPIKE_MUTATION,
        control_iso_sql=WT_VARIANTS,
        pos=pos,
        aa=aa
    )

    _format_print_rx_info(run_query(VP_SQL))


def query_outbreak_variant(pos, aa):
    print('\n==== Outbreak variants ====')
    print('  |  '.join(['Pangolin', 'total', 'URL']))

    variant_count_map = load_variant_count()

    outbreak_info = load_csv(OUTBREAK_PATH)
    for row in outbreak_info:
        if int(row['pos']) == pos and row['mut'].upper() == aa.upper():
            total = variant_count_map.get(row['name'], 'Unknown')
            print('  |  '.join(
                [
                    row['name'],
                    str(total),
                    'https://outbreak.info/situation-reports?pango={}'.format(
                        row['name'])
                ]))
            continue


def _query(pos, aa):
    print('-'*20)
    print("AA:", aa)
    query_dms_by_mutation(pos, aa)
    query_invitro_by_mutation(pos, aa)
    query_invivo_by_mutation(pos, aa)
    query_outbreak_variant(pos, aa)
    print('-'*20 + '\n\n')


def query_muation(pos, aa):
    print('Position:', pos)
    if aa and len(aa) == 1:
        _query(pos, aa)
    else:
        for aa in AA_NAMES:
            _query(pos, aa)


if __name__ == '__main__':
    query_muation(POSITION, AA)
