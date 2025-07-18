from django.shortcuts import render
from math import log2, log10

from process.models import (
    Project,
    SampleStage,
    Abundance,
    Statistic,
    Protein,
    Phospho,
)

from process.constants import (
    ABUNDANCES_RAW,
    ABUNDANCES_NORMALISED_LOG2_MEAN,
    ANOVA,
    Q_VALUE,
    CURVE_FOLD_CHANGE,
    PROJECT_SL,
)

def index(request):
    return render(request, "index.html", {})


def venn(request):
    SL_project = Project.objects.get(name='SL')
    Original_project = Project.objects.get(name='Original')

    SL_proteins = SL_project.phospho_protein_list
    Original_proteins = Original_project.phospho_protein_list

    intersection = set(SL_proteins).intersection(set(Original_proteins))

    sets = [
        {'sets': ['SL'], 'size': len(SL_proteins)},
        {'sets': ['Original'], 'size': len(Original_proteins)},
        {'sets': ['SL', 'Original'], 'size': len(intersection)},
    ]

    return render(request, "venn.html", {
        "venn_data": sets,
    })



def heatmap(request):
    accession_numbers = ['Q8TF01', 'Q8TF05', 'Q8TF42', 'Q8TF68', 'Q8TF72']

    sample_stages = SampleStage.objects.filter(project__name=PROJECT_SL).order_by('rank')

    samples = [s.name.removeprefix("RO sample ").replace("sample, ", "") for s in sample_stages]

    readings = []

    for accession_number in accession_numbers:
        abundances = Abundance.objects.filter(
            statistic__statistic_type__name = ABUNDANCES_RAW,
            statistic__protein__accession_number = accession_number,
            statistic__protein__project__name = PROJECT_SL,
            replicate__mean = True
        ).order_by(
            'sample_stage__rank'
        )

        row = [a.reading for a in abundances]

        readings.append(row)

    return render(request, "heatmap.html", {
        "accession_numbers_data": accession_numbers,
        "samples_data": samples,
        "readings_data": readings,
    })



def pca(request):
    # return render(request, "pca.html", {"pca-data": pca_json})
    return render(request, "pca.html", {})


def barchart(request):
    return render(request, "barchart.html", {})


def scatterplot(request):
    protein_data = []
    phospho_data = []

    project = Project.objects.get(name=PROJECT_SL)

    protein_num = Protein.objects.filter(project=project).count()
    phospho_num = Phospho.objects.filter(protein__project=project).count()

    oscillating_protein_num = len(project.protein_list)
    oscillating_phospho_num = len(project.phospho_protein_list)

    protein_statistics = Statistic.objects.filter(
        statistic_type__name = ABUNDANCES_NORMALISED_LOG2_MEAN,
        protein__project__name = PROJECT_SL,
    )

    phospho_statistics = Statistic.objects.filter(
        statistic_type__name = ABUNDANCES_NORMALISED_LOG2_MEAN,
        phospho__protein__project__name = PROJECT_SL,
    )

    protein_data = process_scatterplot_statistics(protein_statistics)
    phospho_data = process_scatterplot_statistics(phospho_statistics)

    return render(request, "scatterplot.html", {
        "protein_data": protein_data,
        "protein_num": protein_num,
        "oscillating_protein_num": oscillating_protein_num,
        "phospho_data": phospho_data,
        "phospho_num": phospho_num,
        "oscillating_phospho_num": oscillating_phospho_num,
    })

def process_scatterplot_statistics(statistics):
    data = []

    for statistic in statistics:
        if (
            statistic.metrics.get(ANOVA) is not None and
            statistic.metrics[ANOVA].get(Q_VALUE) and
            statistic.metrics.get(CURVE_FOLD_CHANGE)
        ):
            try:
                curve_fold_log2 = log2(abs(statistic.metrics[CURVE_FOLD_CHANGE]))
                q_value_negative_log10 = -log10(statistic.metrics[ANOVA][Q_VALUE])
            except Exception as e:
                print(e)
                print(statistic.metrics[CURVE_FOLD_CHANGE])
                print(statistic.metrics[ANOVA][Q_VALUE])
                continue

            if curve_fold_log2 < 0:
                continue

            # .263 is log2(1.2), -log2(0.01) is 2
            if curve_fold_log2 > .263 and q_value_negative_log10 > 2:
                oscillating = True
            else:
                oscillating = False

            data.append(
                {
                    "x": curve_fold_log2,
                    "y": q_value_negative_log10,
                    "label": "",
                    "oscillating": oscillating,
                }
            )

    return data