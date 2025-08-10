from django.urls import path
from .views import (
    healthz, pdb_proxy, get_protein_info, prepare_ligand, prepare_receptor,
    estimate_binding_site, run_docking, get_docking_status, list_docking_jobs,
    detect_binding_pockets, rescore_poses, run_advanced_docking, 
    list_job_templates, get_pocket_list
)

urlpatterns = [
    path('healthz', healthz, name='healthz'),
    path('pdb/<str:pdb_id>', pdb_proxy, name='pdb-proxy'),
    path('pdb/<str:pdb_id>/info', get_protein_info, name='get-protein-info'),
    path('ligand/prepare', prepare_ligand, name='prepare-ligand'),
    path('receptor/prepare', prepare_receptor, name='prepare-receptor'),
    path('binding-site/estimate', estimate_binding_site, name='estimate-binding-site'),
    path('dock/run', run_docking, name='run-docking'),
    path('dock/status/<str:job_id>', get_docking_status, name='get-docking-status'),
    path('dock/jobs', list_docking_jobs, name='list-docking-jobs'),
    
    # Advanced Features
    path('dock/advanced/run', run_advanced_docking, name='run-advanced-docking'),
    path('pockets/detect', detect_binding_pockets, name='detect-binding-pockets'),
    path('pockets/<str:pdb_id>', get_pocket_list, name='get-pocket-list'),
    path('poses/rescore', rescore_poses, name='rescore-poses'),
    path('templates', list_job_templates, name='list-job-templates'),
]