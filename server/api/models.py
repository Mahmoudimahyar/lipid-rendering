from django.db import models
import uuid
import json
from datetime import datetime
from django.contrib.auth.models import User

class DockingJob(models.Model):
    """Model to track molecular docking jobs"""
    
    # Job identification
    job_id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    
    # Job status
    STATUS_CHOICES = [
        ('pending', 'Pending'),
        ('running', 'Running'),
        ('completed', 'Completed'),
        ('failed', 'Failed'),
        ('cancelled', 'Cancelled'),
    ]
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='pending')
    
    # Input parameters
    ligand_smiles = models.TextField(help_text="SMILES string of the ligand")
    receptor_pdb_id = models.CharField(max_length=10, help_text="PDB ID of the receptor")
    
    # Docking parameters
    center_x = models.FloatField(help_text="X coordinate of binding site center")
    center_y = models.FloatField(help_text="Y coordinate of binding site center")
    center_z = models.FloatField(help_text="Z coordinate of binding site center")
    size_x = models.FloatField(help_text="Search space size in X dimension")
    size_y = models.FloatField(help_text="Search space size in Y dimension")
    size_z = models.FloatField(help_text="Search space size in Z dimension")
    
    exhaustiveness = models.IntegerField(default=8, help_text="Thoroughness of search")
    num_modes = models.IntegerField(default=9, help_text="Number of binding modes to generate")
    seed = models.IntegerField(null=True, blank=True, help_text="Random seed for reproducible results")
    
    # Engine metadata
    engine = models.CharField(max_length=20, null=True, blank=True, help_text="Docking engine used (vina, mock, etc.)")
    is_mock = models.BooleanField(default=False, help_text="Whether mock docking was used")
    
    # Results
    results_json = models.JSONField(null=True, blank=True, help_text="Docking results in JSON format")
    error_message = models.TextField(null=True, blank=True, help_text="Error message if job failed")
    
    # Advanced features
    user_session = models.CharField(max_length=255, null=True, blank=True, help_text="User session identifier")
    job_name = models.CharField(max_length=255, null=True, blank=True, help_text="User-provided job name")
    job_description = models.TextField(null=True, blank=True, help_text="Job description")
    advanced_params = models.JSONField(null=True, blank=True, help_text="Advanced docking parameters")
    rescoring_results = models.JSONField(null=True, blank=True, help_text="GNINA rescoring results")
    pocket_detection_results = models.JSONField(null=True, blank=True, help_text="Binding pocket detection results")
    
    # Performance metrics
    cpu_time = models.FloatField(null=True, blank=True, help_text="CPU time in seconds")
    memory_usage = models.IntegerField(null=True, blank=True, help_text="Peak memory usage in MB")
    # Progress tracking (Phase 4)
    progress_percent = models.FloatField(null=True, blank=True, help_text="Job progress percent (0-100)")
    estimated_seconds_remaining = models.FloatField(null=True, blank=True, help_text="Estimated seconds remaining")
    last_progress_at = models.DateTimeField(null=True, blank=True, help_text="Last progress update timestamp")
    progress_message = models.CharField(max_length=255, null=True, blank=True, help_text="Human-readable progress message")
    
    # Timestamps
    created_at = models.DateTimeField(auto_now_add=True)
    started_at = models.DateTimeField(null=True, blank=True)
    completed_at = models.DateTimeField(null=True, blank=True)
    
    # Tags and categories
    tags = models.JSONField(default=list, help_text="User-defined tags for job organization")
    priority = models.IntegerField(default=1, help_text="Job priority (1=low, 5=high)")
    
    class Meta:
        ordering = ['-created_at']
    
    def __str__(self):
        return f"DockingJob {self.job_id} - {self.status}"
    
    def duration(self):
        """Calculate job duration if completed"""
        if self.started_at and self.completed_at:
            return (self.completed_at - self.started_at).total_seconds()
        return None
    
    def is_finished(self):
        """Check if job is in a terminal state"""
        return self.status in ['completed', 'failed', 'cancelled']    
    def get_best_pose(self):
        """Get the best scoring pose from results"""
        if not self.results_json or 'poses' not in self.results_json:
            return None
        
        poses = self.results_json['poses']
        if not poses:
            return None
            
        # Return pose with best (most negative) affinity
        return min(poses, key=lambda x: x.get('affinity', 0))
    
    def get_tags_display(self):
        """Get comma-separated tags for display"""
        return ', '.join(self.tags) if self.tags else ''


class BindingPocket(models.Model):
    """Model to store detected binding pockets"""
    
    pocket_id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    protein_pdb_id = models.CharField(max_length=10, help_text="PDB ID of the protein")
    
    # Pocket geometry
    center_x = models.FloatField(help_text="X coordinate of pocket center")
    center_y = models.FloatField(help_text="Y coordinate of pocket center") 
    center_z = models.FloatField(help_text="Z coordinate of pocket center")
    radius = models.FloatField(null=True, blank=True, help_text="Pocket radius in Angstroms")
    volume = models.FloatField(help_text="Pocket volume in cubic Angstroms")
    surface_area = models.FloatField(null=True, blank=True, help_text="Pocket surface area in square Angstroms")
    
    # Pocket properties
    druggability_score = models.FloatField(null=True, blank=True, help_text="Druggability score (0-1)")
    hydrophobicity = models.FloatField(null=True, blank=True, help_text="Hydrophobicity score")
    polarity = models.FloatField(null=True, blank=True, help_text="Polarity score")
    
    # Enhanced pocket analysis
    pocket_rank = models.IntegerField(null=True, blank=True, help_text="Rank among detected pockets")
    shape_analysis = models.JSONField(null=True, blank=True, help_text="Shape and geometry analysis")
    chemical_environment = models.JSONField(null=True, blank=True, help_text="Chemical environment analysis")
    druggability_analysis = models.JSONField(null=True, blank=True, help_text="Detailed druggability analysis")
    
    # Detection metadata
    detection_method = models.CharField(max_length=50, default='fpocket', help_text="Method used for detection")
    detection_software = models.CharField(max_length=100, null=True, blank=True, help_text="Software used for detection")
    confidence_score = models.FloatField(help_text="Confidence in pocket detection")
    residues = models.JSONField(help_text="List of residues forming the pocket")
    
    # Additional properties for enhanced analysis
    cavity_type = models.CharField(max_length=50, null=True, blank=True, help_text="Type of cavity (deep, shallow, tunnel, cleft)")
    accessibility = models.CharField(max_length=50, null=True, blank=True, help_text="Accessibility level (buried, surface, partially_buried)")
    conservation_score = models.FloatField(null=True, blank=True, help_text="Evolutionary conservation score")
    
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    
    class Meta:
        ordering = ['-druggability_score', '-confidence_score']
        unique_together = ['protein_pdb_id', 'center_x', 'center_y', 'center_z']
    
    def __str__(self):
        return f"Pocket {self.pocket_id} in {self.protein_pdb_id} (Score: {self.druggability_score})"
    
    def get_center_coordinates(self):
        """Get pocket center as list"""
        return [self.center_x, self.center_y, self.center_z]
    
    def distance_to_point(self, x, y, z):
        """Calculate distance from pocket center to given point"""
        import math
        return math.sqrt((self.center_x - x)**2 + (self.center_y - y)**2 + (self.center_z - z)**2)
    
    def is_druggable(self):
        """Check if pocket is considered druggable"""
        return self.druggability_score and self.druggability_score > 0.5
    
    def get_druggability_class(self):
        """Get druggability classification"""
        if not self.druggability_score:
            return 'unknown'
        elif self.druggability_score > 0.7:
            return 'high'
        elif self.druggability_score > 0.4:
            return 'medium'
        else:
            return 'low'


class JobTemplate(models.Model):
    """Model to store reusable docking job templates"""
    
    template_id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=255, help_text="Template name")
    description = models.TextField(help_text="Template description")
    
    # Template parameters
    default_params = models.JSONField(help_text="Default docking parameters")
    advanced_params = models.JSONField(null=True, blank=True, help_text="Advanced parameters")
    
    # Usage tracking
    usage_count = models.IntegerField(default=0, help_text="Number of times template was used")
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    
    # Categories
    category = models.CharField(max_length=100, default='general', help_text="Template category")
    is_public = models.BooleanField(default=False, help_text="Whether template is publicly available")
    
    class Meta:
        ordering = ['-usage_count', 'name']
    
    def __str__(self):
        return f"Template: {self.name}"
    
    def increment_usage(self):
        """Increment usage counter"""
        self.usage_count += 1
        self.save(update_fields=['usage_count'])
