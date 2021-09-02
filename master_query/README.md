## Run the script

1. synchronize new data set
    - please `cd` to `covid-drdb-reports`
    - please run `make syncdb`
    - `you only need to run this command once in a day.`
2. update the position and mutation
    - open file `master_query_file.py`
    - update `POSITION` and `AA`
3. run `make query` to view the result


## Criterias used in the query

- All control variants are wildtype or B.1 variant
- Vaccine plasmas are from fully vaccinated subjects
- mAb are in EUA or in clinical trials
