RX_FULLY_VACCINE = """
SELECT
    *
FROM
    rx_vacc_plasma
WHERE
    (vaccine_name, dosage) = ('AZD1222', 2)
    OR
    (vaccine_name, dosage) = ('Ad26.COV2.S', 1)
    OR
    (vaccine_name, dosage) = ('Ad26.COV2.S', 2)
    OR
    (vaccine_name, dosage) = ('BBIBP-CorV', 2)
    OR
    (vaccine_name, dosage) = ('BBV152', 2)
    OR
    (vaccine_name, dosage) = ('BBV154', 1)
    OR
    (vaccine_name, dosage) = ('BNT162b2', 2)
    OR
    (vaccine_name, dosage) = ('CoronaVac', 2)
    OR
    (vaccine_name, dosage) = ('CoronaVac', 3)
    OR
    (vaccine_name, dosage) = ('MVC-COV1901', 2)
    OR
    (vaccine_name, dosage) = ('NVX-CoV2373', 2)
    OR
    (vaccine_name, dosage) = ('Sputnik V', 2)
    OR
    (vaccine_name, dosage) = ('ZF2001', 2)
    OR
    (vaccine_name, dosage) = ('ZF2001', 3)
    OR
    (vaccine_name, dosage) = ('mRNA-1273', 2)
    OR
    (vaccine_name, dosage) = ('mRNA-1273', 3)
"""
