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
    projects = Project.objects.all()

    projects_data = []

    for project in projects:
        name = project.name

        if name == "Original":
            name = "ICR from original json output"
        elif name == "ICR":
            if project.with_bugs:
                name = "ICR (with bugs)"
            else:
                name = "ICR (without bugs)"

        projects_data.append({
            "id": project.id,
            "name": name,
        })

    return render(request, "index.html", {
        "projects": projects
    })


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



def heatmap(request, id):
    project = Project.objects.get(id=id)

    # CDK1, CDK2, CDK6, CCND1
    protein_genes = ['CDK1', 'CDK2', 'CDK6', 'CCND1']
    protein_accession_numbers = ['P06493', 'P24941', 'Q00534', 'P24385']

    sample_stages = SampleStage.objects.filter(project=project).order_by('rank')

    samples = [s.name.removeprefix("RO sample ").replace("sample, ", "") for s in sample_stages]

    protein_readings = []
    protein_reading_max = 0

    for accession_number in protein_accession_numbers:
        abundances = Abundance.objects.filter(
            statistic__statistic_type__name = ABUNDANCES_RAW,
            statistic__protein__accession_number = accession_number,
            statistic__protein__project = project,
            replicate__mean = True
        ).order_by(
            'sample_stage__rank'
        )

        row = []

        for a in abundances:
            row.append(a.reading)

            if a.reading > protein_reading_max:
                protein_reading_max = a.reading                

        protein_readings.append(row)


    phospho_phosphosites = ['pY6', 'pT4', 'pT5']
    phospho_peptides = ['IGEGTY(UniMod:21)GVVYK2', 'ALGT(UniMod:21)PNNEVWPEVESLQDYK3', 'IGEGT(UniMod:21)YGVVYK2']

    phosphos = Phospho.objects.filter(
        protein__accession_number = "P06493",
        protein__project = project
    )

    phospho_readings = []
    phospho_reading_max = 0

    for phospho in phosphos:
        abundances = Abundance.objects.filter(
            statistic__statistic_type__name = ABUNDANCES_RAW,
            statistic__phospho = phospho,
            replicate__mean = True
        ).order_by(
            'sample_stage__rank'
        )

        row = []

        for a in abundances:
            row.append(a.reading)

            if a.reading > phospho_reading_max:
                phospho_reading_max = a.reading                

        phospho_readings.append(row)

    return render(request, "heatmap.html", {
        "protein_genes": protein_genes,
        "protein_readings_data": protein_readings,
        "protein_reading_max": protein_reading_max,
        "phospho_phosphosites": [p.mod for p in phosphos],
        "phospho_readings_data": phospho_readings,
        "phospho_reading_max": phospho_reading_max,
        "samples_data": samples,
    })



def pca(request):
    # return render(request, "pca.html", {"pca-data": pca_json})
    return render(request, "pca.html", {})


def barchart(request, id):
    project = Project.objects.get(id = id)

    protein_go = project.protein_go_list
    phospho_go = project.phospho_protein_go_list

    protein_bar_data = process_bar_data(protein_go)
    phospho_bar_data = process_bar_data(phospho_go)

    return render(request, "barchart.html", {
        "protein_data": protein_bar_data,
        "phospho_data": phospho_bar_data,
    })

def process_bar_data(go_list):
    top = {
        'CC': [],
        'BP': [],
        'MF': []
    }

    groups = {
        'GO:CC': 'CC',
        'GO:BP': 'BP',
        'GO:MF': 'MF',
    }

    if go_list:
        for go in go_list:
            source = go['source']
            p_value = go['p_value']

            if len(top[groups[source]]) < 3:
                top[groups[source]].append({
                    "group": groups[source],
                    "name": go["name"],
                    "value": -log10(p_value)
                })

    return [item for sublist in top.values() for item in sublist]


def scatterplot(request, id):
    project = Project.objects.get(id=id)

    protein_data = []
    phospho_data = []

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