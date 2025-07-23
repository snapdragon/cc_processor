### For local development
```sh
eval $(poetry env activate)
poetry install
```


### To run from scratch
```sh
docker compose up
docker exec -it cc_processor /bin/bash
# poetry install        ?
eval $(poetry env activate)
python manage.py import_protein --project "SL"
python manage.py import_phospho --project "SL"
python manage.py process --project "SL" --calculate-all

python manage.py import_protein --project "ICR"
python manage.py import_phospho --project "ICR"
python manage.py process --project "ICR" --with-bugs --calculate-all

# Or the long route
python manage.py process --project "ICR" --with-bugs --calculate-protein-medians
python manage.py process --project "ICR" --with-bugs --calculate-proteins
python manage.py process --project "ICR" --with-bugs --calculate-phospho-medians
python manage.py process --project "ICR" --with-bugs --calculate-phosphos
python manage.py process --project "ICR" --with-bugs --calculate-batch
python manage.py process --project "ICR" --with-bugs --fetch-references

```

### To run with IRC
```sh
python manage.py import_spreadsheet --project "ICR"
python manage.py process --project "ICR"
```

### To reset the migration
```sh
docker exec -it postgres-db /bin/bash
psql -U myuser -d mydatabase
\dt
# DROP TABLE 'process_' for all tables starting 'process', e.g.
drop table process_peptidestartposition;
drop table process_uniprotdata;
drop table process_abundance;
drop table process_statistic;
drop table process_statistictype;
drop table process_phospho;
drop table process_protein;
drop table process_columnname;
drop table process_samplestage;
drop table process_replicate;
drop table process_project;
delete from django_migrations where app = 'process';
# Exit the
# Go to the code (in vscode or other editor)
# Copy the contents of the migration file with the inserts, e.g. from 0002_auto_20250531_1657.
# Delete the migration files
# Go to the process docker container
docker exec -it cc_processor /bin/bash
python manage.py makemigrations
python manage.py makemigrations --empty process
# Paste all the previously copied inserts into the file made by the command above
python manage.py migrate
# From inside postgres container
psql -U myuser -d dbify -f db_data/process_peptidestartposition.sql
psql -U myuser -d dbify -f db_data/process_uniprotdata.sql

```

### Run docker compose
```sh
docker compose up
```

### Connect to django container
```sh
docker exec -it cc_processor /bin/bash
```

### Create and run migrations
```sh
python manage.py makemigrations
python manage.py migrate
python manage.py makemigrations --empty process
```

### Connect to postgres container
```sh
docker exec -it postgres-db /bin/bash
psql -U myuser -d mydatabase
```

### How to clear out a migration
```sh
delete from django_migrations where app = 'process';
```
### Run venv
```sh
eval $(poetry env activate)
```

### Build docker image
```sh
docker build -t cc_processor .
```

### Run the container
```sh
docker run -it --rm cc_processor
```

### Run all pre-commit hooks
```sh
pre-commit run --all-files
```


### Running tests
```sh
pytest
pytest --cov=process
```

### Get json result from DB (example)
```sh
docker exec -it postgres-db /bin/bash
\copy (select protein_result from process_runresult, process_protein where protein_id = process_protein.id and accession_number = 'Q09666' and run_id = 2) TO 'Q09666_postgres.json';
# Now exit container and copy file to local using the command below
docker cp postgres-db:/Q09666_postgres.json ./
```

### Output the results for a protein to the /output dir

```sh
python manage.py output_protein --project "ICR" --with-bugs --accession-number Q93075
```

### Various DB queries
```sh
CREATE DATABASE dbify;

# Change DB
\c dbify

# List tables
\dt

# Display the length of non-empty phospho results in order
SELECT id, protein_id, LENGTH(phospho_result::text) AS json_text_length FROM process_runresult where phospho_result != '{}' order by json_text_length asc;

```

### Import the original ICR output into the DB
```sh
python manage.py import_original
```

### Compare a protein from the original output with new output, or compare all
```sh
python manage.py compare --accession-number Q93075
python manage.py compare
```

### To get the admin site working
```sh
docker exec -it cc_processor /bin/bash
python manage.py createsuperuser
python manage.py runserver 0.0.0.0:8000
# Then go to 0.0.0.0:8000/admin in your browser
```

### Dump tables to save retrieving again
```sh
# From inside postgres_db container
pg_dump -d main -t process_uniprotdata > uniprot_data.sql
pg_dump -d main -t process_peptidestartposition > peptidestartposition.sql

# From outside container
docker cp f09f37276ef6:/process_uniprotdata.sql db_data/process_uniprotdata.sql
docker cp f09f37276ef6:/process_peptidestartposition.sql db_data/process_peptidestartposition.

# To import, from outside container
psql -U myuser -d dbify -f db_data/process_peptidestartposition.sql
psql -U myuser -d dbify -f db_data/process_uniprotdata.sql
```
