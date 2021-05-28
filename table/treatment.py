EXPERIMENTS_VIEW = """
    susc_results
"""


SAMPLES_VIEW = """
(
    SELECT
        DISTINCT ref_name, rx_name, ordinal_number, cumulative_count
    FROM
        susc_results
)
"""
