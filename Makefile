tables:
	@python table/gen_table_data.py ../covid-drdb/local/covid-drdb-latest.db

syncdb:
	@curl https://api.github.com/repos/hivdb/covid-drdb-payload/releases/latest | grep "browser_download_url" | cut -d : -f 2,3 | tr -d \" | sort | tail -1  | wget -i - -O master_query/covid-drdb-latest.db
	@wget https://raw.githubusercontent.com/hivdb/chiro-cms/master/resources/outbreak-aapcnt/variants-mutations.csv -O master_query/variants-mutations.csv

query:
	@python master_query/master_query_file.py > master_query/query_results.txt
	@echo "Please open \033[0;31mmaster_query/query_results.txt\033[0m to review the result."

.PHONE: tables syncdb query
