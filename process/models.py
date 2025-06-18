from django.db import models

# Create your models here.


class Project(models.Model):
    name = models.CharField(max_length=255)
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

    def __str__(self):
        return self.name


class SampleStage(models.Model):
    name = models.CharField(max_length=255)
    rank = models.IntegerField()
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
        return self.sample_stage.name


class Protein(models.Model):
    project = models.ForeignKey(
        Project, on_delete=models.CASCADE, related_name="protein"
    )
    accession_number = models.CharField(max_length=255)
    is_contaminant = models.BooleanField(default=False)

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

    def __str__(self):
        return f"{self.reading}"


class Phospho(models.Model):
    protein = models.ForeignKey(
        Protein, on_delete=models.CASCADE, related_name="phospho"
    )

    # TODO - is this long enough?
    mod = models.CharField(max_length=255)

    # TODO - what is this for?
    phosphosite = models.CharField(max_length=255)

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

    def __str__(self):
        return f"{self.reading}"


class Run(models.Model):
    project = models.ForeignKey(
        Project, on_delete=models.CASCADE, related_name="run"
    )
    limit_proteins = models.BooleanField(default=False)
    with_bugs = models.BooleanField(default=False)
    protein_medians = models.JSONField(blank=True, null=True)
    phospho_medians = models.JSONField(blank=True, null=True)
    results = models.JSONField(blank=True, null=True)

    def __str__(self):
        return f"{self.project.name} {self.limit_proteins}"

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
    protein_phospho_result = models.JSONField(blank=True, null=True)

