## CC_processor

Code to process proteomic and phosphoproteomic cell cycle stage data.
To run, the relevant files must be placed in the /data directory.

# To add new datasets

Create a new Project record for the project, and new Replicate, SampleStage and ColumnNames records for the names of the replicates, stages and columns in the spreadsheet(s). For examples, see the migration files.
You may need to alter import_protein and import_phospho depending on the type of spreadsheet.


### To run from scratch
```sh
# Build cc_processor
docker build -t cc_processor .

docker compose up

# Go to the cc_processor container and set things up
docker exec -it cc_processor /bin/bash
eval $(poetry env activate)
python manage.py migrate

# From outside the container, import the uniprot and peptide data tables
# (not necessary but saves time)
docker exec -i postgres-db psql -U myuser -d main < db_data/process_uniprotdata_2025_07_19.sql
docker exec -i postgres-db psql -U myuser -d main < db_data/process_peptidestartposition_2025_07_19.sql

# From inside the cc_processor container, import and process the data
python manage.py import_protein --project "SL"
python manage.py import_phospho --project "SL"
python manage.py process --project "SL" --run-all
```

### Connect to django container
```sh
docker exec -it cc_processor /bin/bash
```

### Connect to postgres container
```sh
docker exec -it postgres-db /bin/bash
psql -U myuser -d main
```

### To run the webserver (useful for viewing charts)
```sh
docker exec -it cc_processor /bin/bash
python manage.py runserver 0.0.0.0:8000
# In a browser go to http://0.0.0.0:8000/
```

### Import the original ICR output into the DB, useful for comparison
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

### Dump entire DB
``` sh
pg_dump -U myuser -d main \
  -t process_project -t process_replicate -t process_samplestage \
  -t process_columnname -t process_protein -t process_phospho \
  -t process_statistictype -t process_statistic -t process_abundance -t \
  -t process_uniprotdata -t process_peptidestartposition \
  --data-only -F c -f db_backup.sql

# From outside container
docker cp postgres-db:/db_backup.sql db_backups/db_backup.sql

# To import again, first go to the postgres container and log in to psql
# Then run these statements. You will see various errors, ignore them
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

# Run the django migration from inside cc_processor
python manage.py migrate

# Then exit the postgres-db container and run this locally
docker exec -i postgres-db psql -U myuser -d main < db_backups/db_backup.sql
```

### Dump tables to save retrieving again
```sh
# From inside postgres_db container
pg_dump -d main -t process_uniprotdata > uniprot_data.sql
pg_dump -d main -t process_peptidestartposition > peptidestartposition.sql

# From outside container
docker cp postgres-db:/process_uniprotdata.sql db_data/process_uniprotdata.sql
docker cp postgres-db:/process_peptidestartposition.sql db_data/process_peptidestartposition.

# From outside container
docker cp db_data/process_uniprotdata.sql postgres-db:/process_uniprotdata.sql 
docker cp db_data/process_peptidestartposition.sql postgres-db:/process_peptidestartposition.sql 

# To import, from outside container
psql -U myuser -d dbify -f db_data/process_peptidestartposition.sql
psql -U myuser -d dbify -f db_data/process_uniprotdata.sql
```

### Various DB queries
```sh
# Display the length of non-empty phospho results in order
SELECT id, protein_id, LENGTH(phospho_result::text) AS json_text_length FROM process_runresult where phospho_result != '{}' order by json_text_length asc;

# Count phosphos for a project
SELECT p.id AS project_id,
       p.name AS project_name,
       COUNT(ph.id) AS phospho_count
FROM process_project p
JOIN process_protein pr ON pr.project_id = p.id
JOIN process_phospho ph ON ph.protein_id = pr.id
GROUP BY p.id, p.name;

```
