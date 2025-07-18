from django.contrib import admin
from django.urls import path
from process import views

urlpatterns = [
    path('', views.index, name='index'),
    path('venn/', views.venn, name='venn'),
    path('heatmap/<int:id>/', views.heatmap, name='heatmap'),
    path('pca/<int:id>/', views.pca, name='pca'),
    path('barchart/<int:id>/', views.barchart, name='barchart'),
    path('scatterplot/<int:id>/', views.scatterplot, name='scatterplot'),
    path("admin/", admin.site.urls),
]
