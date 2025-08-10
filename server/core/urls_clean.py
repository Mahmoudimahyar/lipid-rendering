"""
URL configuration for core project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include
from django.http import JsonResponse
from django.shortcuts import redirect

def root_view(request):
    """Root endpoint providing API information"""
    return JsonResponse({
        'name': 'Lipid Rendering API',
        'version': '1.0.0',
        'description': 'Molecular docking and visualization platform',
        'endpoints': {
            'health': '/api/healthz',
            'protein_info': '/api/pdb/{pdb_id}/info',
            'docking': '/api/dock/run',
            'advanced_docking': '/api/dock/advanced/run',
            'pocket_detection': '/api/pockets/detect',
            'templates': '/api/templates',
            'admin': '/admin/',
            'frontend': 'http://localhost:3000'
        },
        'status': 'running',
        'features': [
            'Basic molecular docking',
            'Advanced GNINA rescoring',
            'Binding pocket detection',
            'Job templates',
            'Interactive visualization',
            'Multi-pose analysis'
        ]
    })

urlpatterns = [
    path('', root_view, name='root'),
    path('admin/', admin.site.urls),
    path('api/', include('api.urls')),
]
