from django.contrib import admin
from django.urls import path

from process import views

urlpatterns = [
    path("", views.index, name="index"),
    path("venn/", views.venn, name="venn"),
    path("heatmap/<int:id>/", views.heatmap, name="heatmap"),
    path("heatmap/", views.heatmap, name="heatmap"),
    path(
        "heatmap_by_protein/<int:id>/",
        views.heatmap_by_protein,
        name="heatmap_by_protein",
    ),
    path("pca/<int:id>/", views.pca, name="pca"),
    path("barchart/<int:id>/", views.barchart, name="barchart"),
    path("scatterplot/<int:id>/", views.scatterplot, name="scatterplot"),
    path("mean_gene_effect/<int:id>/", views.mean_gene_effect, name="mean_gene_effect"),
    path("half_life/<int:id>/", views.half_life, name="half_life"),
    path("lingress/<int:id>/", views.lingress, name="lingress"),
    path("admin/", admin.site.urls),
]
