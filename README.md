### Run docker compose
```sh
docker compose up
```

### Connect to container
```sh
docker exec -it cc_processor /bin/bash
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

### Run venv
```sh
eval $(poetry env activate)
poetry install
```
