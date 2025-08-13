"""
CI-specific Django settings for Lipid Rendering server

This configuration is optimized for CI environments where:
- Scientific libraries may not be available initially
- Mock mode should be allowed for testing
- Testing should be fast and reliable
"""

from .settings import *
import sys

# Production settings - NO MOCK ALLOWED
DOCKING_ALLOW_MOCK = False  # NEVER use mock in any environment
DOCKING_FORCE_REAL = True   # Always require real AutoDock Vina
DOCKING_CUDA_ENABLED = True  # Enable CUDA when available

# Use file-based database for persistence (required for Docker)
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': '/app/db_ci.sqlite3',
    }
}

# Disable migrations during testing for speed
class DisableMigrations:
    def __contains__(self, item):
        return True
    def __getitem__(self, item):
        return None

if 'test' in sys.argv or 'pytest' in sys.modules:
    MIGRATION_MODULES = DisableMigrations()

# Logging configuration for CI

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'verbose': {
            'format': '{levelname} {asctime} {module} {process:d} {thread:d} {message}',
            'style': '{',
        },
        'simple': {
            'format': '{levelname} {message}',
            'style': '{',
        },
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'simple',
            'stream': sys.stdout,
        },
        'file': {
            'class': 'logging.FileHandler',
            'formatter': 'verbose',
            'filename': 'ci_test.log',
        },
    },
    'root': {
        'handlers': ['console', 'file'],
        'level': 'INFO',
    },
    'loggers': {
        'django': {
            'handlers': ['console'],
            'level': 'WARNING',
            'propagate': False,
        },
        'api': {
            'handlers': ['console', 'file'],
            'level': 'INFO',
            'propagate': False,
        },
    },
}

# Security settings for CI (less strict for testing)
SECRET_KEY = 'ci-test-key-not-for-production'
DEBUG = True
ALLOWED_HOSTS = ['localhost', '127.0.0.1', '0.0.0.0']

# CORS settings for CI testing
CORS_ALLOW_ALL_ORIGINS = True
CORS_ALLOW_CREDENTIALS = True

# Disable CSRF for CI API testing
CSRF_TRUSTED_ORIGINS = [
    'http://localhost:8000',
    'http://127.0.0.1:8000',
    'http://0.0.0.0:8000',
]

# Cache configuration (use dummy cache for CI)
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.dummy.DummyCache',
    }
}

# Static files (simplified for CI)
STATIC_ROOT = '/app/staticfiles'

# Email backend (use console for CI)
EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'

# Time zone for CI
USE_TZ = True
TIME_ZONE = 'UTC'
