from django.http import JsonResponse, HttpResponseBadRequest, HttpResponse
from django.views.decorators.http import require_GET, require_POST
from django.views.decorators.csrf import csrf_exempt
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import status
from django.utils import timezone
import json
import logging
import requests
import threading

from .models import DockingJob, BindingPocket, JobTemplate
from .chem_utils import ChemUtils
from .docking_utils import DockingEngine
from .advanced_docking import AdvancedDockingEngine, PocketDetector, GNINAScorer, JobTemplateManager

logger = logging.getLogger(__name__)

@require_GET
def healthz(request):
    """Health check endpoint"""
    return JsonResponse({"ok": True})

@require_GET
def pdb_proxy(request, pdb_id: str):
    """Proxy endpoint to fetch PDB structure"""
    if not pdb_id or len(pdb_id) != 4:
        return HttpResponseBadRequest("Invalid PDB ID. Must be 4 characters.")
    
    try:
        # Get protein information
        protein_info = ChemUtils.get_protein_info(pdb_id)
        if not protein_info:
            return HttpResponseBadRequest(f"PDB {pdb_id} not found")
        
        # Fetch the actual structure
        pdb_content = ChemUtils.fetch_protein_structure(pdb_id)
        if not pdb_content:
            return HttpResponseBadRequest(f"Failed to fetch PDB structure for {pdb_id}")
        
        # Estimate binding sites
        binding_sites = ChemUtils.estimate_binding_site(pdb_content)
        
        return JsonResponse({
            "pdb_id": pdb_id.upper(),
            "title": protein_info.get('title', ''),
            "organism": protein_info.get('organism', []),
            "method": protein_info.get('method', ''),
            "resolution": protein_info.get('resolution', []),
            "molecular_weight": protein_info.get('molecular_weight', 0),
            "deposited": protein_info.get('deposited', ''),
            "pdb_content": pdb_content,
            "binding_sites": binding_sites.get('binding_sites', []),
            "fetch_method": "RCSB PDB API"
        })
    
    except Exception as e:
        logger.error(f"PDB proxy error: {e}")
        return HttpResponseBadRequest(f"Internal server error: {e}")

@require_GET
def get_protein_info(request, pdb_id: str):
    """Get protein metadata from RCSB PDB"""
    try:
        info = ChemUtils.get_protein_info(pdb_id)
        return JsonResponse(info or {})
    except Exception as e:
        return HttpResponseBadRequest(f"Failed to fetch protein info: {e}")

@csrf_exempt
@require_POST
def prepare_ligand(request):
    """Prepare ligand for docking from SMILES input"""
    try:
        body = json.loads(request.body or "{}")
    except json.JSONDecodeError:
        return HttpResponseBadRequest("Invalid JSON")

    smiles = (body.get("smiles") or "").strip()
    num_conformers = int(body.get("num_conformers", 1))
    if not smiles:
        return HttpResponseBadRequest("SMILES is required")

    try:
        result = ChemUtils.prepare_ligand_for_docking(smiles)
        return JsonResponse(result)
    except ValueError as e:
        return HttpResponseBadRequest(str(e))
    except requests.exceptions.RequestException as e:
        return HttpResponseBadRequest(f"Failed to prepare ligand: {e}")

@csrf_exempt
@require_POST
def prepare_receptor(request):
    """Prepare receptor for docking from PDB input"""
    try:
        body = json.loads(request.body or "{}")
    except json.JSONDecodeError:
        return HttpResponseBadRequest("Invalid JSON")

    pdb_id = (body.get("pdb_id") or "").strip()
    if not pdb_id:
        return HttpResponseBadRequest("PDB ID is required")

    try:
        result = ChemUtils.prepare_receptor_for_docking(pdb_id)
        return JsonResponse(result)
    except ValueError as e:
        return HttpResponseBadRequest(str(e))
    except requests.exceptions.RequestException as e:
        return HttpResponseBadRequest(f"Failed to prepare receptor: {e}")

@csrf_exempt
@require_POST
def estimate_binding_site(request):
    """Estimate binding site for a given PDB structure"""
    try:
        body = json.loads(request.body or "{}")
    except json.JSONDecodeError:
        return HttpResponseBadRequest("Invalid JSON")

    pdb_id = (body.get("pdb_id") or "").strip()
    ligand_smiles = (body.get("ligand_smiles") or "").strip()
    
    if not pdb_id:
        return HttpResponseBadRequest("PDB ID is required")

    try:
        binding_site = DockingEngine.estimate_binding_site_auto(pdb_id, ligand_smiles or None)
        return JsonResponse(binding_site)
    except Exception as e:
        logger.error(f"Binding site estimation error: {e}")
        return HttpResponseBadRequest(f"Failed to estimate binding site: {e}")

@csrf_exempt
@require_POST
def run_docking(request):
    """
    Start a molecular docking job
    """
    try:
        body = json.loads(request.body or "{}")
    except json.JSONDecodeError:
        return HttpResponseBadRequest("Invalid JSON")

    try:
        # Validate docking parameters
        validated_params = DockingEngine.validate_docking_parameters(body)
        
        # Create docking job
        job = DockingJob.objects.create(
            ligand_smiles=validated_params['ligand_smiles'],
            receptor_pdb_id=validated_params['receptor_pdb_id'],
            center_x=validated_params['center_x'],
            center_y=validated_params['center_y'],
            center_z=validated_params['center_z'],
            size_x=validated_params['size_x'],
            size_y=validated_params['size_y'],
            size_z=validated_params['size_z'],
            exhaustiveness=validated_params['exhaustiveness'],
            num_modes=validated_params['num_modes'],
            status='pending'
        )
        
        # Start docking in background
        threading.Thread(target=_run_docking_task, args=(job.job_id,), daemon=True).start()
        
        return JsonResponse({
            'job_id': str(job.job_id),
            'status': job.status,
            'message': 'Docking job started',
            'estimated_time': '2-5 minutes'
        })
        
    except ValueError as e:
        return HttpResponseBadRequest(str(e))
    except Exception as e:
        logger.error(f"Docking job creation error: {e}")
        return HttpResponseBadRequest(f"Failed to start docking job: {e}")

@require_GET
def get_docking_status(request, job_id: str):
    """Get status of a docking job"""
    try:
        job = DockingJob.objects.get(job_id=job_id)
        
        response_data = {
            'job_id': str(job.job_id),
            'status': job.status,
            'created_at': job.created_at.isoformat(),
            'ligand_smiles': job.ligand_smiles,
            'receptor_pdb_id': job.receptor_pdb_id,
        }
        
        if job.started_at:
            response_data['started_at'] = job.started_at.isoformat()
        
        if job.completed_at:
            response_data['completed_at'] = job.completed_at.isoformat()
            response_data['duration'] = job.duration()
        
        if job.error_message:
            response_data['error_message'] = job.error_message
            
        if job.results_json:
            response_data['results'] = job.results_json
            
        return JsonResponse(response_data)
        
    except DockingJob.DoesNotExist:
        return HttpResponseBadRequest(f"Docking job {job_id} not found")
    except Exception as e:
        logger.error(f"Get docking status error: {e}")
        return HttpResponseBadRequest(f"Failed to get job status: {e}")

@require_GET
def list_docking_jobs(request):
    """List recent docking jobs"""
    try:
        # Get last 20 jobs
        jobs = DockingJob.objects.all()[:20]
        
        jobs_data = []
        for job in jobs:
            job_data = {
                'job_id': str(job.job_id),
                'status': job.status,
                'ligand_smiles': job.ligand_smiles,
                'receptor_pdb_id': job.receptor_pdb_id,
                'created_at': job.created_at.isoformat(),
            }
            
            if job.completed_at:
                job_data['completed_at'] = job.completed_at.isoformat()
                job_data['duration'] = job.duration()
                
            jobs_data.append(job_data)
        
        return JsonResponse({
            'jobs': jobs_data,
            'total_count': len(jobs_data)
        })
        
    except Exception as e:
        logger.error(f"List docking jobs error: {e}")
        return HttpResponseBadRequest(f"Failed to list jobs: {e}")

def _run_docking_task(job_id):
    """
    Background task to run docking calculation
    """
    try:
        job = DockingJob.objects.get(job_id=job_id)
        
        # Update job status to running
        job.status = 'running'
        job.started_at = timezone.now()
        job.save()
        
        # Prepare docking parameters
        validated_params = {
            'ligand_smiles': job.ligand_smiles,
            'receptor_pdb_id': job.receptor_pdb_id,
            'center_x': job.center_x,
            'center_y': job.center_y,
            'center_z': job.center_z,
            'size_x': job.size_x,
            'size_y': job.size_y,
            'size_z': job.size_z,
            'exhaustiveness': job.exhaustiveness,
            'num_modes': job.num_modes,
        }
        
        # Prepare inputs
        inputs = DockingEngine.prepare_docking_inputs(
            job.ligand_smiles,
            job.receptor_pdb_id
        )
        
        # Run docking calculation (real or mock)
        results = DockingEngine.run_production_docking(validated_params)
        
        # Update job with results
        job.status = 'completed'
        job.completed_at = timezone.now()
        job.results_json = results
        job.save()
        
        logger.info(f"Docking job {job_id} completed successfully")
        
    except Exception as e:
        # Update job with error
        try:
            job = DockingJob.objects.get(job_id=job_id)
            job.status = 'failed'
            job.completed_at = timezone.now()
            job.error_message = str(e)
            job.save()
        except:
            pass
        
        logger.error(f"Docking job {job_id} failed: {e}")


# Advanced Features Endpoints

@csrf_exempt
@require_POST
def detect_binding_pockets(request):
    """
    Detect binding pockets in a protein structure
    """
    try:
        body = json.loads(request.body or "{}")
    except json.JSONDecodeError:
        return HttpResponseBadRequest("Invalid JSON")

    pdb_id = (body.get("pdb_id") or "").strip()
    if not pdb_id:
        return HttpResponseBadRequest("PDB ID is required")

    try:
        # Get protein structure
        pdb_content = ChemUtils.fetch_protein_structure(pdb_id)
        if not pdb_content:
            return HttpResponseBadRequest(f"Could not fetch structure for PDB {pdb_id}")

        # Detect pockets
        detected_pockets = PocketDetector.detect_pockets(pdb_id, pdb_content)
        
        # Store pockets in database
        saved_pockets = []
        for pocket_data in detected_pockets:
            pocket, created = BindingPocket.objects.get_or_create(
                protein_pdb_id=pdb_id,
                center_x=pocket_data['center_x'],
                center_y=pocket_data['center_y'],
                center_z=pocket_data['center_z'],
                defaults={
                    'radius': pocket_data['radius'],
                    'volume': pocket_data['volume'],
                    'druggability_score': pocket_data['druggability_score'],
                    'hydrophobicity': pocket_data['hydrophobicity'],
                    'polarity': pocket_data['polarity'],
                    'detection_method': pocket_data['detection_method'],
                    'confidence_score': pocket_data['confidence_score'],
                    'residues': pocket_data['residues']
                }
            )
            
            # Add druggability analysis
            pocket_data['druggability_analysis'] = PocketDetector.analyze_pocket_druggability(pocket_data)
            saved_pockets.append(pocket_data)

        return JsonResponse({
            'pdb_id': pdb_id,
            'pockets': saved_pockets,
            'summary': {
                'total_pockets': len(saved_pockets),
                'druggable_pockets': len([p for p in saved_pockets if p['druggability_score'] > 0.6]),
                'detection_method': 'fpocket'
            }
        })
        
    except Exception as e:
        logger.error(f"Pocket detection error: {e}")
        return HttpResponseBadRequest(f"Failed to detect pockets: {e}")


@csrf_exempt
@require_POST
def rescore_poses(request):
    """
    Rescore docking poses using GNINA
    """
    try:
        body = json.loads(request.body or "{}")
    except json.JSONDecodeError:
        return HttpResponseBadRequest("Invalid JSON")

    job_id = (body.get("job_id") or "").strip()
    if not job_id:
        return HttpResponseBadRequest("Job ID is required")

    try:
        job = DockingJob.objects.get(job_id=job_id)
        
        if not job.results_json or 'poses' not in job.results_json:
            return HttpResponseBadRequest("Job has no poses to rescore")

        poses = job.results_json['poses']
        
        # Rescore poses using GNINA
        rescored_poses = GNINAScorer.rescore_poses(poses, "mock_sdf", "mock_pdbqt")
        
        # Update job with rescoring results
        job.rescoring_results = {
            'method': 'GNINA',
            'rescored_poses': rescored_poses,
            'timestamp': timezone.now().isoformat()
        }
        job.save()

        return JsonResponse({
            'job_id': str(job.job_id),
            'rescoring_method': 'GNINA',
            'rescored_poses': rescored_poses,
            'summary': {
                'original_best_score': min(poses, key=lambda x: x['affinity'])['affinity'],
                'gnina_best_score': min(rescored_poses, key=lambda x: x['gnina_score'])['gnina_score'],
                'poses_rescored': len(rescored_poses)
            }
        })
        
    except DockingJob.DoesNotExist:
        return HttpResponseBadRequest(f"Docking job {job_id} not found")
    except Exception as e:
        logger.error(f"Pose rescoring error: {e}")
        return HttpResponseBadRequest(f"Failed to rescore poses: {e}")


@csrf_exempt
@require_POST
def run_advanced_docking(request):
    """
    Start an advanced molecular docking job with additional features
    """
    try:
        body = json.loads(request.body or "{}")
    except json.JSONDecodeError:
        return HttpResponseBadRequest("Invalid JSON")

    try:
        # Validate basic docking parameters
        validated_params = DockingEngine.validate_docking_parameters(body)
        
        # Extract advanced options
        use_pocket_detection = body.get('use_pocket_detection', False)
        use_gnina_rescoring = body.get('use_gnina_rescoring', False)
        job_name = (body.get('job_name') or '').strip()
        job_description = (body.get('job_description') or '').strip()
        tags = body.get('tags', [])
        template_id = (body.get('template_id') or '').strip()
        
        # Apply template if specified
        if template_id:
            try:
                template = JobTemplate.objects.get(template_id=template_id)
                template_params = JobTemplateManager.apply_template(template, validated_params)
                validated_params.update(template_params)
            except JobTemplate.DoesNotExist:
                return HttpResponseBadRequest(f"Template {template_id} not found")
        
        # Create advanced docking job
        job = DockingJob.objects.create(
            ligand_smiles=validated_params['ligand_smiles'],
            receptor_pdb_id=validated_params['receptor_pdb_id'],
            center_x=validated_params['center_x'],
            center_y=validated_params['center_y'],
            center_z=validated_params['center_z'],
            size_x=validated_params['size_x'],
            size_y=validated_params['size_y'],
            size_z=validated_params['size_z'],
            exhaustiveness=validated_params['exhaustiveness'],
            num_modes=validated_params['num_modes'],
            job_name=job_name or f"Advanced Docking {validated_params['receptor_pdb_id']}",
            job_description=job_description,
            tags=tags,
            advanced_params={
                'use_pocket_detection': use_pocket_detection,
                'use_gnina_rescoring': use_gnina_rescoring,
                'template_id': template_id
            },
            status='pending'
        )
        
        # Start advanced docking in background
        threading.Thread(
            target=_run_advanced_docking_task, 
            args=(job.job_id, use_pocket_detection, use_gnina_rescoring), 
            daemon=True
        ).start()
        
        return JsonResponse({
            'job_id': str(job.job_id),
            'status': job.status,
            'job_name': job.job_name,
            'message': 'Advanced docking job started',
            'features': {
                'pocket_detection': use_pocket_detection,
                'gnina_rescoring': use_gnina_rescoring
            },
            'estimated_time': '3-8 minutes'
        })
        
    except ValueError as e:
        return HttpResponseBadRequest(str(e))
    except Exception as e:
        logger.error(f"Advanced docking job creation error: {e}")
        return HttpResponseBadRequest(f"Failed to start advanced docking job: {e}")


@require_GET
def list_job_templates(request):
    """
    List available job templates
    """
    try:
        category = request.GET.get('category', '')
        
        if category:
            templates = JobTemplate.objects.filter(category=category, is_public=True)
        else:
            templates = JobTemplate.objects.filter(is_public=True)
        
        templates_data = []
        for template in templates:
            templates_data.append({
                'template_id': str(template.template_id),
                'name': template.name,
                'description': template.description,
                'category': template.category,
                'usage_count': template.usage_count,
                'default_params': template.default_params,
                'advanced_params': template.advanced_params
            })
        
        return JsonResponse({
            'templates': templates_data,
            'categories': list(JobTemplate.objects.values_list('category', flat=True).distinct())
        })
        
    except Exception as e:
        logger.error(f"List templates error: {e}")
        return HttpResponseBadRequest(f"Failed to list templates: {e}")


@require_GET
def get_pocket_list(request, pdb_id: str):
    """
    Get stored binding pockets for a protein
    """
    try:
        pockets = BindingPocket.objects.filter(protein_pdb_id=pdb_id.upper())
        
        pockets_data = []
        for pocket in pockets:
            pocket_data = {
                'pocket_id': str(pocket.pocket_id),
                'center_x': pocket.center_x,
                'center_y': pocket.center_y,
                'center_z': pocket.center_z,
                'radius': pocket.radius,
                'volume': pocket.volume,
                'druggability_score': pocket.druggability_score,
                'confidence_score': pocket.confidence_score,
                'detection_method': pocket.detection_method,
                'created_at': pocket.created_at.isoformat()
            }
            pockets_data.append(pocket_data)
        
        return JsonResponse({
            'pdb_id': pdb_id.upper(),
            'pockets': pockets_data,
            'total_pockets': len(pockets_data)
        })
        
    except Exception as e:
        logger.error(f"Get pocket list error: {e}")
        return HttpResponseBadRequest(f"Failed to get pocket list: {e}")


def _run_advanced_docking_task(job_id, use_pocket_detection=False, use_gnina_rescoring=False):
    """
    Background task to run advanced docking calculation
    """
    try:
        job = DockingJob.objects.get(job_id=job_id)
        
        # Update job status to running
        job.status = 'running'
        job.started_at = timezone.now()
        job.save()
        
        # Prepare docking parameters
        validated_params = {
            'ligand_smiles': job.ligand_smiles,
            'receptor_pdb_id': job.receptor_pdb_id,
            'center_x': job.center_x,
            'center_y': job.center_y,
            'center_z': job.center_z,
            'size_x': job.size_x,
            'size_y': job.size_y,
            'size_z': job.size_z,
            'exhaustiveness': job.exhaustiveness,
            'num_modes': job.num_modes,
        }
        
        # Run advanced docking calculation
        results = AdvancedDockingEngine.run_advanced_docking(
            validated_params,
            use_pocket_detection=use_pocket_detection,
            use_gnina_rescoring=use_gnina_rescoring
        )
        
        # Update job with results and performance metrics
        job.status = 'completed'
        job.completed_at = timezone.now()
        job.results_json = results
        
        if 'performance' in results:
            job.cpu_time = results['performance'].get('cpu_time')
            job.memory_usage = results['performance'].get('memory_usage')
        
        if 'detected_pockets' in results:
            job.pocket_detection_results = results['detected_pockets']
        
        if 'rescored_poses' in results:
            job.rescoring_results = {
                'method': 'GNINA',
                'rescored_poses': results['rescored_poses']
            }
        
        job.save()
        
        logger.info(f"Advanced docking job {job_id} completed successfully")
        
    except Exception as e:
        # Update job with error
        try:
            job = DockingJob.objects.get(job_id=job_id)
            job.status = 'failed'
            job.completed_at = timezone.now()
            job.error_message = str(e)
            job.save()
        except:
            pass
        
        logger.error(f"Advanced docking job {job_id} failed: {e}")
