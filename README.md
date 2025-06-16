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
python manage.py import_proteo --project "Soliman Labs"
python manage.py import_phospho --project "Soliman Labs"
python manage.py process --project "Soliman Labs"

python manage.py import_proteo --project "ICR"
python manage.py import_phospho --project "ICR"
python manage.py process --project "ICR" --with-bugs
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
drop table process_phosphoreading;
drop table process_phospho;
drop table process_proteinreading;
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