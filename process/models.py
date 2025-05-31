from django.db import models

# Create your models here.


class Project(models.Model):
    name = models.CharField(max_length=255)
    proteome_file = models.CharField(max_length=255)
    phosphoproteome_file = models.CharField(max_length=255)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.name


class Replicate(models.Model):
    name = models.CharField(max_length=255)
    rank = models.IntegerField
    project = models.ForeignKey(
        Project, on_delete=models.CASCADE, related_name="replicates"
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.name


class SampleStage(models.Model):
    name = models.CharField(max_length=255)
    rank = models.IntegerField
    project = models.ForeignKey(
        Project, on_delete=models.CASCADE, related_name="sample_stages"
    )

    def __str__(self):
        return self.name


class ColumnName(models.Model):
    name = models.CharField(max_length=255)
    sample_stage = models.ForeignKey(
        SampleStage, on_delete=models.CASCADE, related_name="column_names"
    )
    replicate = models.ForeignKey(
        Replicate, on_delete=models.CASCADE, related_name="sample_stages"
    )

    def __str__(self):
        return self.name
