# Production Dockerfile for Helix Guard Backend
FROM python:3.12-slim

WORKDIR /app

# Install system dependencies needed for scientific libraries & PDF generation
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy and install dependencies
COPY requirements.txt requirements.txt
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy backend files
COPY . .

# Set environment
ENV PYTHONUNBUFFERED=1
ENV PORT=8000
ENV PYTHONPATH=/app

EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
  CMD curl -f http://127.0.0.1:${PORT}/api/health || exit 1

# Launch uvicorn server
CMD ["sh", "-c", "uvicorn backend.app.main:app --host 0.0.0.0 --port ${PORT} || uvicorn app.main:app --host 0.0.0.0 --port ${PORT}"]
