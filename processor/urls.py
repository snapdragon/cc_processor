from django.contrib import admin
from django.urls import path
from process import views

urlpatterns = [
    path('', views.index, name='index'),
    path('venn/', views.venn, name='venn'),
    path('heatmap/', views.heatmap, name='heatmap'),
    path('pca/', views.pca, name='pca'),
    path('barchart/', views.barchart, name='barchart'),
    path('scatterplot/', views.scatterplot, name='scatterplot'),
    path("admin/", admin.site.urls),
]
