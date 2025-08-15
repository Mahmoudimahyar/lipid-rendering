# Lipid Rendering Server

This directory contains the Django backend for the Lipid Rendering platform. It provides a robust REST API for molecular docking, data management, and real-time visualization support.

## 🧬 Features

- **Django REST Framework**: For building a powerful and flexible API.
- **AutoDock Vina Integration**: For performing real molecular docking.
- **Celery and Redis**: For asynchronous task management and handling long-running docking jobs.
- **Scientific Libraries**: Utilizes RDKit, OpenBabel, and Meeko for cheminformatics tasks.
- **Advanced Docking Features**: Includes support for binding pocket detection and pose rescoring with GNINA.
- **Dockerized Environment**: Fully containerized for consistent development and production environments.

## 🏗️ Architecture

The backend is structured into two main Django apps:

- **`core`**: The main project app, responsible for settings, URL routing, and serving the React frontend.
- **`api`**: The core application that contains all the business logic for molecular docking, including models, views, and utility functions.

## 🔧 Environment Configuration

Copy the following into a `.env` file (or set as environment variables):

```
DJANGO_SECRET_KEY=change-me
DEBUG=False
ALLOWED_HOSTS=*
CORS_ALLOWED_ORIGINS=http://localhost:3000,http://127.0.0.1:3000
CSRF_TRUSTED_ORIGINS=http://localhost:3000,http://127.0.0.1:3000
LOG_LEVEL=INFO
WHITENOISE_MAX_AGE=31536000
```

### Key Components

- **`api/views.py`**: Contains the view functions that handle incoming API requests and orchestrate the docking process.
- **`api/models.py`**: Defines the Django models for storing docking jobs, results, and other relevant data.
- **`api/docking_utils.py`**: A utility module that encapsulates the core logic for running docking jobs and processing results.
- **`api/real_docking_engine.py`**: Manages the integration with AutoDock Vina and other scientific libraries.
- **`core/urls.py`**: The main URL configuration file that routes requests to the appropriate views.

## 🚀 Getting Started

### Prerequisites

- **Docker and Docker Compose**: For running the application in a containerized environment.
- **Python 3.11+**: For local development and running management commands.

### Local Development

1.  **Navigate to the project root**:
    ```bash
    cd lipid-rendering
    ```

2.  **Build and start the Docker containers**:
    ```bash
    docker-compose up --build
    ```

3.  **The backend server will be running at `http://localhost:8000`**.

### Running Tests

To run the backend tests, you can execute the following command from the project root:

```bash
docker-compose run --rm server pytest
```

## 🌐 API Endpoints

The backend exposes a comprehensive set of API endpoints for various functionalities. For a detailed list of all available endpoints, please refer to the `api/urls.py` file or the main project `README.md`.

- **Health Check**: `GET /api/healthz`
- **Docking Capabilities**: `GET /api/dock/capabilities`
- **Run Docking**: `POST /api/dock/run`
- **Get Docking Status**: `GET /api/dock/status/<job_id>`

## 🤝 Contributing

Please refer to the main `CONTRIBUTING.md` file in the project root for guidelines on how to contribute to the backend development.
