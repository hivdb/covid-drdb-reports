tables:
	@python table/gen_table_data.py ../covid-drdb/local/covid-drdb-latest.db

.PHONE: tables
