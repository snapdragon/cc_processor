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
        return self.name


class Protein(models.Model):
    accession_number = models.CharField(max_length=255)

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


# class ProcessResult(models.Model):
#     column_name = models.ForeignKey(
#         ColumnName, on_delete=models.CASCADE, related_name="protein_readings"
#     )
#     created_at = models.DateTimeField(auto_now_add=True)
#     updated_at = models.DateTimeField(auto_now=True)

#     def __str__(self):
#         return self.column_name.sample_stage.name
