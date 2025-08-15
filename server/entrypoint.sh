#!/usr/bin/env bash
set -euo pipefail

# Run migrations quietly; ignore if no DB changes
python manage.py migrate --noinput || true

# Collect static (idempotent)
python manage.py collectstatic --noinput || true

exec gunicorn core.wsgi:application --bind 0.0.0.0:8000 --workers 2


