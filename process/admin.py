from django.contrib import admin

from .models import Abundance, Statistic, StatisticType, Phospho, Protein, ColumnName, SampleStage, Replicate, Project

admin.site.register(Abundance)
admin.site.register(Statistic)
admin.site.register(StatisticType)
admin.site.register(Phospho)
admin.site.register(Protein)
admin.site.register(ColumnName)
admin.site.register(SampleStage)
admin.site.register(Replicate)
admin.site.register(Project)
