### To run from scratch
```sh
docker compose up
docker exec -it cc_processor /bin/bash
eval $(poetry env activate)
python processor/scripts/import_spreadsheet.py
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
