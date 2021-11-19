from statistics import stdev
from statistics import mean


def get_z_score(num_list):
    if len(num_list) == 0:
        raise Exception('Z-Score: the list is null')

    mean_value = mean(num_list)
    std_value = stdev(num_list)

    result = []
    for num in num_list:
        z_score = (num - mean_value) / std_value
        result.append((num, z_score))

    return result


def get_outlier(record_list, key):

    if len(record_list) < 2:
        return []

    num_list = [i[key] for i in record_list]
    mean_value = mean(num_list)
    std_value = stdev(num_list)

    if not std_value:
        return []

    result = []
    for rec in record_list:
        num = rec[key]
        z_score = (num - mean_value) / std_value
        if z_score > 3:
            result.append(rec)

    return result

