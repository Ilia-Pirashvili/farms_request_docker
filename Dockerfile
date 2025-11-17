FROM python:3.11-slim

WORKDIR /app

# Install system dependencies for rasterio, Pillow, etc.
RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential \
  gdal-bin \
  libgdal-dev \
  libpq-dev \
  libexpat1 \
  && rm -rf /var/lib/apt/lists/*

# Copy Python dependencies first (for caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

ENV HDF5_USE_FILE_LOCKING=FALSE

# Copy app code
COPY app/ ./app

# Copy all files (background_data, user_data, etc.)
COPY files/ ./files

# Copy authenticators
RUN mkdir -p /root/.local/share/openeo-python-client/
## Copernicus
COPY files/persistant/refresh-tokens.json /root/.local/share/openeo-python-client/refresh-tokens.json

# Expose volume for optional override
VOLUME /files

# Run the FastAPI app
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", "--timeout-keep-alive", "0"]
