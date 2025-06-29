from django.db import models

# Create your models here.


class Project(models.Model):
    name = models.CharField(max_length=255, unique=True)
    proteome_file = models.CharField(max_length=255)
    phosphoproteome_file = models.CharField(max_length=255)
    proteome_file_accession_number_column_name = models.CharField(max_length=255)

    def __str__(self):
        return self.name

class Replicate(models.Model):
    name = models.CharField(max_length=255)
    rank = models.IntegerField()
    project = models.ForeignKey(
        Project, on_delete=models.CASCADE, related_name="replicates"
    )

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['project', 'name'], name='unique_replicate_project_name')
        ]

    def __str__(self):
        return self.name

class SampleStage(models.Model):
    name = models.CharField(max_length=255)
    rank = models.IntegerField()
    project = models.ForeignKey(
        Project, on_delete=models.CASCADE, related_name="sample_stages"
    )

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['project', 'name'], name='unique_sample_stage_project_name')
        ]

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

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['name', 'replicate', 'sample_stage'], name='unique_name_replicate_sample_stage')
        ]

    def __str__(self):
        return self.sample_stage.name


class Protein(models.Model):
    project = models.ForeignKey(
        Project, on_delete=models.CASCADE, related_name="protein"
    )
    accession_number = models.CharField(max_length=255)
    is_contaminant = models.BooleanField(default=False)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['project', 'accession_number'], name='unique_project_accession_number')
        ]

    def __str__(self):
        return self.accession_number


class ProteinReading(models.Model):
    column_name = models.ForeignKey(
        ColumnName, on_delete=models.CASCADE, related_name="protein_readings"
    )
    protein = models.ForeignKey(
        Protein, on_delete=models.CASCADE, related_name="protein_readings"
    )
    reading = models.FloatField(null=True, blank=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['column_name', 'protein'], name='unique_column_name_protein')
        ]

    def __str__(self):
        return f"{self.reading}"


class Phospho(models.Model):
    protein = models.ForeignKey(
        Protein, on_delete=models.CASCADE, related_name="phospho"
    )

    mod = models.CharField(max_length=255)

    phosphosite = models.CharField(max_length=255)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['protein', 'mod'], name='unique_protein_mod')
        ]

    def __str__(self):
        return f"Phospho for {self.protein.accession_number}"


class PhosphoReading(models.Model):
    column_name = models.ForeignKey(
        ColumnName, on_delete=models.CASCADE, related_name="phospho_reading"
    )
    phospho = models.ForeignKey(
        Phospho, on_delete=models.CASCADE, related_name="phospho_reading"
    )
    reading = models.FloatField(null=True, blank=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['phospho', 'column_name'], name='unique_phospho_column_name')
        ]

    def __str__(self):
        return f"{self.reading}"


class Run(models.Model):
    project = models.ForeignKey(
        Project, on_delete=models.CASCADE, related_name="run"
    )
    with_bugs = models.BooleanField(default=False)
    protein_medians = models.JSONField(blank=True, null=True)
    phospho_medians = models.JSONField(blank=True, null=True)
    # TODO - is this field used?
    results = models.JSONField(blank=True, null=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['project', 'with_bugs'], name='unique_project_with_bugs')
        ]

    def __str__(self):
        return f"{self.project.name}, with-bugs {self.with_bugs}"

# TODO - duplication here - proteins are per-project, so both run and protein link to protein
#   Does it matter? Not sure if it's fixable.
class RunResult(models.Model):
    run = models.ForeignKey(
        Run, on_delete=models.CASCADE, related_name="run_result"
    )
    protein = models.ForeignKey(
        Protein, on_delete=models.CASCADE, related_name="run_result"
    )
    protein_result = models.JSONField(blank=True, null=True)
    phospho_result = models.JSONField(blank=True, null=True)
    combined_result = models.JSONField(blank=True, null=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['run', 'protein'], name='unique_run_protein')
        ]
