FROM python:3.12-slim

# Don't write .pyc files
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

ENV POETRY_VERSION=2.0.1

WORKDIR /app

# Set environment variables to avoid prompts and caching issues
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Update system packages and install pip dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip to the latest version
RUN python -m ensurepip && \
    pip install --upgrade pip

# Set working directory
WORKDIR /app

# Default command
CMD ["python3"]
