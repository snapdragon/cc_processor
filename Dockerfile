FROM python:3.12-slim

# Environment settings
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    POETRY_VERSION=2.0.1 \
    PATH="/app/.venv/bin:$PATH"

WORKDIR /app

# Install system dependencies, including R and R headers
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl build-essential gcc gfortran \
    libpq-dev libreadline-dev \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libbz2-dev zlib1g-dev liblzma-dev \
    ca-certificates \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Poetry
RUN pip install "poetry==$POETRY_VERSION"

# Copy and install Python dependencies
COPY pyproject.toml poetry.lock ./
RUN poetry config virtualenvs.in-project true && \
    poetry install --only main --no-root

# Copy application source code
COPY . .

# Expose Django development port
EXPOSE 8000

# Set the default command
CMD ["/bin/bash"]
