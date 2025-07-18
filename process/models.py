from django.db import models
from django.db.models import Q, UniqueConstraint

class Project(models.Model):
    name = models.CharField(max_length=255, unique=True)

    with_bugs = models.BooleanField(default=False)
    # Whether the project has abundance data that can be processed.
    #   Currently only applies to ICR and SL.
    processable = models.BooleanField(default=False)
    proteome_file = models.CharField(max_length=255)
    phosphoproteome_file = models.CharField(max_length=255)
    proteome_file_accession_number_column_name = models.CharField(max_length=255)
    protein_list = models.JSONField(null=True, blank=True)
    phospho_protein_list = models.JSONField(null=True, blank=True)
    # GO is lowercase as postgres automatically folds unquoted identifiers to lowercase
    protein_go_list =  models.JSONField(null=True, blank=True)
    phospho_protein_go_list =  models.JSONField(null=True, blank=True)

    def __str__(self):
        return f"Project {self.name}"

class Replicate(models.Model):
    name = models.CharField(max_length=255)
    mean = models.BooleanField(default=False)
    project = models.ForeignKey(
        Project, on_delete=models.CASCADE, related_name="replicate_project"
    )

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['project', 'name'], name='unique_replicate_project_name')
        ]

    def __str__(self):
        return f"Replicate {self.name}"

class SampleStage(models.Model):
    name = models.CharField(max_length=255)
    rank = models.IntegerField()
    project = models.ForeignKey(
        Project, on_delete=models.CASCADE, related_name="sample_stage_project"
    )

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['project', 'name'], name='unique_sample_stage_project_name')
        ]

    def __str__(self):
        return f"SampleStage {self.name}"


class ColumnName(models.Model):
    name = models.CharField(max_length=255)
    sample_stage = models.ForeignKey(
        SampleStage, on_delete=models.CASCADE, related_name="column_name_sample_stage"
    )
    replicate = models.ForeignKey(
        Replicate, on_delete=models.CASCADE, related_name="sample_stage_replicate"
    )

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['name', 'replicate', 'sample_stage'], name='unique_name_replicate_sample_stage')
        ]

    def __str__(self):
        return f"ColumnName {self.sample_stage.name}"


class Protein(models.Model):
    project = models.ForeignKey(
        Project, on_delete=models.CASCADE, related_name="protein_project"
    )
    accession_number = models.CharField(max_length=255)
    is_contaminant = models.BooleanField(default=False)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['project', 'accession_number'], name='unique_project_accession_number')
        ]

    def __str__(self):
        return f"Protein {self.accession_number}"


class Phospho(models.Model):
    protein = models.ForeignKey(
        Protein, on_delete=models.CASCADE, related_name="phospho_protein"
    )

    mod = models.CharField(max_length=255)

    phosphosite = models.CharField(max_length=255)

    pep_tools_annotations = models.JSONField(null=True, blank=True)
    kinase_prediction = models.JSONField(null=True, blank=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['protein', 'mod'], name='unique_protein_mod')
        ]

    def __str__(self):
        return f"Phospho for {self.protein.accession_number} mod {self.mod}"


class StatisticType(models.Model):
    name = models.CharField(max_length=255)

    def __str__(self):
        return f"StatisticType {self.name}"


class Statistic(models.Model):
    statistic_type = models.ForeignKey(
        StatisticType, on_delete=models.CASCADE, related_name="statistic_type"
    )
    protein = models.ForeignKey(
        Protein, on_delete=models.CASCADE, related_name="statistic_protein",
        null=True, blank=True
    )
    phospho = models.ForeignKey(
        Phospho, on_delete=models.CASCADE, related_name="statistic_phospho",
        null=True, blank=True
    )
    project = models.ForeignKey(  # Needed for medians
        Project, on_delete=models.SET_NULL, related_name="statistic_project",
        null=True, blank=True
    )
    metrics = models.JSONField(null=True, blank=True)

    class Meta:
        constraints = [
            UniqueConstraint(
                fields=["statistic_type", "project"],
                condition=Q(protein__isnull=True) & Q(phospho__isnull=True),
                name="unique_statistic_by_type_and_project"
            ),
            UniqueConstraint(
                fields=["statistic_type", "protein"],
                condition=Q(project__isnull=True) & Q(phospho__isnull=True),
                name="unique_statistic_by_type_and_protein"
            ),
            UniqueConstraint(
                fields=["statistic_type", "phospho"],
                condition=Q(project__isnull=True) & Q(protein__isnull=True),
                name="unique_statistic_by_type_and_phospho"
            )
        ]

        
    def __str__(self):
        if self.phospho:
            return f"Statistic {self.statistic_type.name} {self.phospho.mod} (phospho)"
        elif self.protein:
            return f"Statistic {self.statistic_type.name} {self.protein.accession_number} (protein)"
        elif self.project:
            return f"Statistic {self.statistic_type.name} {self.project.name} (project)"
        else:
            return f"Invalid statistic {self.statistic_type.name} (no protein or phospho)"


class Abundance(models.Model):
    statistic = models.ForeignKey(
        Statistic, on_delete=models.CASCADE, related_name="abundance_statistic"
    )
    sample_stage = models.ForeignKey(
        SampleStage, on_delete=models.CASCADE, related_name="abundance_sample_stage"
    )
    replicate = models.ForeignKey(
        Replicate, on_delete=models.CASCADE, related_name="abundance_replicate"
    )
    reading = models.FloatField(null=True, blank=True)

    def __str__(self):
        return f"Abundance: replicate '{self.replicate.name}' stage '{self.sample_stage.name}' statistic type '{self.statistic.statistic_type.name}' reading: {self.reading}"

class UniprotData(models.Model):
    accession_number = models.CharField(max_length=255)
    data = models.JSONField(null=True, blank=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['accession_number'], name='unique_uniprot_data_accession_number')
        ]

    def __str__(self):
        return f"Uniprot data for  {self.accession_number}"
