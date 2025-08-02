from math import log2, log10, ceil
import pandas as pd
from sklearn.decomposition import PCA
import json
import numpy as np
from sklearn.preprocessing import StandardScaler

from django.shortcuts import render

from process.constants import (
    ABUNDANCES_NORMALISED_LOG2_MEAN,
    ABUNDANCES_RAW,
    ANOVA,
    CURVE_FOLD_CHANGE,
    PROJECT_SL,
    Q_VALUE,
    ABUNDANCES_NORMALISED_MIN_MAX,
)
from process.models import (
    Abundance,
    Phospho,
    Project,
    Protein,
    SampleStage,
    Statistic,
    UniprotData,
    Replicate,
)


def index(request):
    projects = Project.objects.all()

    return render(request, "index.html", {"projects": projects})


def venn(request):
    SL_project = Project.objects.get(name="SL")
    ICR_project = Project.objects.get(name="ICR")

    SL_proteins = Protein.objects.filter(
        project=SL_project,
        ccd=True
    )
    ICR_proteins = Protein.objects.filter(
        project=ICR_project,
        ccd=True
    )

    SL_accession_numbers = [p.accession_number for p in SL_proteins]
    ICR_accession_numbers = [p.accession_number for p in ICR_proteins]

    intersection = set(SL_accession_numbers).intersection(set(ICR_accession_numbers))

    SL_set_name = f"{len(SL_accession_numbers)}"
    ICR_set_name = f"{len(ICR_accession_numbers)}"

    print(f"intersection {len(intersection)}")

    sets = [
        {"sets": [SL_set_name], "size": len(SL_proteins)},
        {"sets": [ICR_set_name], "size": len(ICR_proteins)},
        {"sets": [SL_set_name, ICR_set_name], "size": len(intersection)},
    ]

    return render(
        request,
        "venn.html",
        {
            "venn_data": sets,
        },
    )


def heatmap_by_protein(request, id):
    project = Project.objects.get(id=id)

    proteins = Protein.objects.filter(project=project)

    return render(
        request,
        "heatmap_by_protein.html",
        {
            "proteins": proteins,
            "project_id": id,
        },
    )


def heatmap(request, id=None):
    if id is None:
        # It's a call from heatmap_by_protein
        project = Project.objects.get(id=request.POST.get("project_id"))

        accession_number = request.POST.get("accession_number")

        protein_accession_numbers = [accession_number]

        try:
            updata = UniprotData.objects.get(accession_number=accession_number)

            protein_genes = [updata.gene_name]
        except Exception:
            # No uniprot data for this accession number
            protein_genes = [accession_number]
    else:
        project = Project.objects.get(id=id)

        protein_genes = ["CDK1", "CDK2", "CDK6", "CCND1"]
        protein_accession_numbers = ["P06493", "P24941", "Q00534", "P24385"]

        accession_number = None

    sample_stages = SampleStage.objects.filter(project=project).order_by("rank")

    samples = [
        s.name.removeprefix("RO sample ").replace("sample, ", "") for s in sample_stages
    ]

    protein_readings = []
    protein_reading_max = 0

    for an in protein_accession_numbers:
        previous_rank = 0

        abundances = Abundance.objects.filter(
            statistic__statistic_type__name=ABUNDANCES_RAW,
            statistic__protein__accession_number=an,
            statistic__protein__project=project,
            replicate__mean=True,
        ).order_by("sample_stage__rank")

        row = []

        for a in abundances:
            while a.sample_stage.rank > previous_rank + 1:
                row.append(None)
                previous_rank += 1

            row.append(a.reading)

            previous_rank = a.sample_stage.rank

            if a.reading > protein_reading_max:
                protein_reading_max = a.reading

        protein_readings.append(row)

    if accession_number is not None:
        phosphos = Phospho.objects.filter(
            protein__accession_number=accession_number, protein__project=project
        )
    else:
        # accession number for CDK1
        phosphos = Phospho.objects.filter(
            protein__accession_number="P06493", protein__project=project
        )

    phospho_readings = []
    phospho_reading_max = 0

    for phospho in phosphos:
        previous_rank = 0

        abundances = Abundance.objects.filter(
            statistic__statistic_type__name=ABUNDANCES_RAW,
            statistic__phospho=phospho,
            replicate__mean=True,
        ).order_by("sample_stage__rank")

        row = []

        for a in abundances:
            while a.sample_stage.rank > previous_rank + 1:
                row.append(None)
                previous_rank += 1

            row.append(a.reading)

            previous_rank = a.sample_stage.rank

            if a.reading > phospho_reading_max:
                phospho_reading_max = a.reading

        phospho_readings.append(row)

    return render(
        request,
        "heatmap.html",
        {
            "project": project,
            "protein_genes": protein_genes,
            "protein_readings_data": protein_readings,
            "protein_reading_max": protein_reading_max,
            "phospho_phosphosites": [p.phosphosite for p in phosphos],
            "phospho_readings_data": phospho_readings,
            "phospho_reading_max": phospho_reading_max,
            "samples_data": samples,
            "accession_number": accession_number,
        },
    )


def pca(request, id):
    project = Project.objects.get(id=id)
    
    return render(request, "pca.html", {
        "protein_data": project.protein_pca,
        "phospho_data": project.phospho_pca
    })


def barchart(request, id):
    project = Project.objects.get(id=id)

    protein_go = project.protein_go_list
    phospho_go = project.phospho_protein_go_list

    protein_bar_data = process_bar_data(protein_go)
    phospho_bar_data = process_bar_data(phospho_go)

    return render(
        request,
        "barchart.html",
        {
            "protein_data": protein_bar_data,
            "phospho_data": phospho_bar_data,
            "project": project,
        },
    )


def process_bar_data(go_list):
    top = {"CC": [], "BP": [], "MF": []}

    groups = {
        "GO:CC": "CC",
        "GO:BP": "BP",
        "GO:MF": "MF",
    }

    if go_list:
        for go in go_list:
            source = go["source"]
            p_value = go["p_value"]

            if len(top[groups[source]]) < 3:
                top[groups[source]].append(
                    {
                        "group": groups[source],
                        "name": go["name"],
                        "value": -log10(p_value),
                    }
                )

    return [item for sublist in top.values() for item in sublist]


def scatterplot(request, id):
    project = Project.objects.get(id=id)

    protein_data = []
    phospho_data = []

    protein_num = (
        Abundance.objects.filter(
            statistic__statistic_type__name=ABUNDANCES_RAW,
            statistic__protein__project=project,
        )
        .values("statistic__protein")
        .distinct()
        .count()
    )

    phospho_num = (
        Abundance.objects.filter(
            statistic__statistic_type__name=ABUNDANCES_RAW,
            statistic__phospho__protein__project=project,
        )
        .values("statistic__phospho")
        .distinct()
        .count()
    )

    oscillating_protein_num = Protein.objects.filter(
        project=project,
        ccd=True
    ).count()
    oscillating_phospho_num = Phospho.objects.filter(
        protein__project=project,
        ccd=True
    ).count()

    protein_statistics = Statistic.objects.filter(
        statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
        protein__project = project,
    )

    phospho_statistics = Statistic.objects.filter(
        statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
        phospho__protein__project = project
    )

    protein_data = process_scatterplot_statistics(protein_statistics, True)
    phospho_data = process_scatterplot_statistics(phospho_statistics, False)

    return render(
        request,
        "scatterplot.html",
        {
            "project": project,
            "protein_data": protein_data,
            "protein_num": protein_num,
            "oscillating_protein_num": oscillating_protein_num,
            "phospho_data": phospho_data,
            "phospho_num": phospho_num,
            "oscillating_phospho_num": oscillating_phospho_num,
        },
    )



def lingress(request, id):
    project = Project.objects.get(id=id)

    sample_stages = SampleStage.objects.filter(project=project).order_by("rank")

    samples = [
        s.name.removeprefix("RO sample ").replace("sample, ", "") for s in sample_stages
    ]

    protein_readings = []
    protein_reading_max = 0
    protein_genes = []

    proteins = Protein.objects.filter(
        project=project,
        ccd=True
    )

    for protein in proteins:
        previous_rank = 0

        protein_genes.append(protein.accession_number)

        abundances = Abundance.objects.filter(
            statistic__statistic_type__name=ABUNDANCES_NORMALISED_MIN_MAX,
            statistic__protein=protein,
            replicate__mean=True,
        ).order_by("sample_stage__rank")

        row = []

        for a in abundances:
            while a.sample_stage.rank > previous_rank + 1:
                row.append(None)
                previous_rank += 1

            row.append(a.reading)

            previous_rank = a.sample_stage.rank

            if a.reading > protein_reading_max:
                protein_reading_max = a.reading

        protein_readings.append(row)

    return render(
        request,
        "lingress.html",
        {
            "project": project,
            "protein_genes": protein_genes,
            "protein_readings_data": protein_readings,
            "protein_reading_max": protein_reading_max,
            "samples_data": samples,
        },
    )





def process_scatterplot_statistics(statistics, is_protein):
    data = []

    important_accession_numbers = ["P24385", "O96020", "P14635", "P20248", "P30281", "P06493", "P24941", "Q00534"]

    for statistic in statistics:
        if (
            statistic.metrics.get(ANOVA) is not None
            and statistic.metrics[ANOVA].get(Q_VALUE)
            and statistic.metrics.get(CURVE_FOLD_CHANGE)
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

            # if is_protein:
            #     # .263 is log2(1.2), -log10(0.01) is 2
            #     if curve_fold_log2 > 0.263 and q_value_negative_log10 > 2:
            #         oscillating = True
            #     else:
            #         oscillating = False
            # else:
            #     # .263 is log2(1.2), -log10(0.01) is 2
            #     if curve_fold_log2 > 0.263 and q_value_negative_log10 > 2:
            #         oscillating = True
            #     else:
            #         oscillating = False

            if is_protein:
                oscillating = statistic.protein.ccd
            else:
                oscillating = statistic.phospho.ccd
                
            gene_name = None

            if is_protein:
                accession_number = statistic.protein.accession_number
            else:
                accession_number = statistic.phospho.protein.accession_number

            if accession_number in important_accession_numbers:
                try:
                    updata = UniprotData.objects.get(accession_number=accession_number)

                    if is_protein:
                        gene_name = updata.gene_name
                    # CCND3-S264 not known to be significant, can't include automatically
                    elif statistic.phospho.phosphosite != "S264":
                        gene_name = f"{updata.gene_name}-{statistic.phospho.phosphosite}"
                except Exception:
                    continue

            data.append(
                {
                    "x": curve_fold_log2,
                    "y": q_value_negative_log10,
                    "label": gene_name,
                    "oscillating": oscillating,
                }
            )

    return data
