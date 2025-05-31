FROM python:3.12-slim

# Don't write .pyc files
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

ENV POETRY_VERSION=2.0.1

WORKDIR /app

# Set environment variables to avoid prompts and caching issues
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
curl build-essential libpq-dev && \
pip install "poetry==$POETRY_VERSION" && \
apt-get clean && rm -rf /var/lib/apt/lists/*

# Upgrade pip to the latest version
RUN python -m ensurepip && \
    pip install --upgrade pip

# Copy only dependency files first
COPY pyproject.toml poetry.lock ./

# Install dependencies into a virtual environment in /opt/venv
RUN poetry config virtualenvs.in-project true && \
    poetry install --only main --no-root

# Copy application code
COPY . .

# Default command
CMD ["/bin/bash"]
