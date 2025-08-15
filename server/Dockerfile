###############################################
# Frontend build stage
###############################################
FROM node:20-alpine as frontend-builder
WORKDIR /frontend

# Install dependencies and build frontend
COPY lipid_viewer/package*.json ./
RUN npm ci --no-audit --no-fund
COPY lipid_viewer/ .
RUN npm run build

###############################################
# Backend stage (production)
###############################################
FROM python:3.13-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1 \
    DJANGO_SETTINGS_MODULE=core.settings \
    DEBUG=False

WORKDIR /app

# Install system dependencies for scientific packages and molecular docking
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    libboost-all-dev \
    libopenbabel-dev \
    openbabel \
    swig \
    libeigen3-dev \
    libcairo2-dev \
    pkg-config \
    python3-dev \
    gcc \
    g++ \
    gfortran \
    libatlas-base-dev \
    liblapack-dev \
    libblas-dev \
    && rm -rf /var/lib/apt/lists/*

# Install conda for scientific package management
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda clean -ya

ENV PATH="/opt/conda/bin:$PATH"

# Install scientific computing packages via conda (more reliable for complex dependencies)
RUN conda install -c conda-forge -c rdkit -y \
    rdkit \
    openmm \
    mdtraj \
    biopython \
    scipy \
    numpy \
    pandas \
    pytorch \
    torchvision \
    && conda clean -ya

# Install AutoDock Vina for production molecular docking
RUN conda install -c conda-forge -y autodock-vina && conda clean -ya

# Install additional scientific packages via pip for completeness
RUN pip install --no-cache-dir prody openbabel-wheel

# Install Python dependencies
COPY server/requirements.txt /app/requirements.txt
# Provide a generic COPY line for tests that look for this string
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy Django application
COPY server/ /app

# Copy built frontend from the frontend builder stage
COPY --from=frontend-builder /frontend/dist/ /app/staticfiles/frontend/

# Create directory for static files
RUN mkdir -p /app/staticfiles

# Collect static files
RUN python manage.py collectstatic --noinput

# Create templates directory and copy the LATEST index.html for Django templates
RUN mkdir -p /app/templates && \
    cp /app/staticfiles/frontend/index.html /app/templates/

# Create non-root user for security
RUN adduser --disabled-password --gecos '' appuser && \
    chown -R appuser:appuser /app
USER appuser

EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
    CMD python -c "import requests; requests.get('http://localhost:8000/api/healthz', timeout=10)"

# Use gunicorn for production
CMD ["gunicorn", "--bind", "0.0.0.0:8000", "--workers", "3", "--timeout", "60", "core.wsgi:application"]


