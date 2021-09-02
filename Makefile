tables:
	@python table/gen_table_data.py ../covid-drdb/local/covid-drdb-latest.db

query:
	@python master_query/master_query_file.py

.PHONE: tables
