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
from django.urls import path, include, re_path
from django.http import JsonResponse, HttpResponse
from django.shortcuts import redirect
from django.views.generic import TemplateView
from django.conf import settings
from django.conf.urls.static import static
from django.views.decorators.cache import never_cache
from django.utils.decorators import method_decorator
import os

def api_info_view(request):
    """API endpoint providing API information"""
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
            'frontend': '/' # Now served from the same server
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

@method_decorator(never_cache, name='dispatch')
class ReactAppView(TemplateView):
    """
    Serve React app for all non-API routes.
    This handles React Router client-side routing.
    """
    template_name = 'index.html'
    
    def get(self, request, *args, **kwargs):
        """Serve React app with no-cache headers to prevent stale HTML"""
        response = super().get(request, *args, **kwargs)
        # Add cache-busting headers for the main HTML
        response['Cache-Control'] = 'no-cache, no-store, must-revalidate'
        response['Pragma'] = 'no-cache'
        response['Expires'] = '0'
        return response
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['api_base_url'] = '/api'
        return context

urlpatterns = [
    # API routes
    path('api/info', api_info_view, name='api-info'),
    path('admin/', admin.site.urls),
    path('api/', include('api.urls')),
]

# Serve static files (including assets at /assets/ path for React compatibility)
# These must come BEFORE the catch-all React route
urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
urlpatterns += static('/assets/', document_root=os.path.join(str(settings.STATIC_ROOT), 'frontend', 'assets'))

# React app routes (catch-all for SPA routing)
# This must be last to catch all non-API/static routes
urlpatterns += [re_path(r'^.*$', ReactAppView.as_view(), name='react-app')]
