tables:
	@python table/work.py ./databases/covid-drdb-latest.db

syncdb:
	@curl https://api.github.com/repos/hivdb/covid-drdb-payload/releases/latest | grep "browser_download_url" | cut -d : -f 2,3 | tr -d \" | sort | tail -1  | wget -i - -O databases/covid-drdb-latest.db

ref_names:
	@python table/ref_names.py ./databases/covid-drdb-latest.db

.PHONE: tables ref_names
