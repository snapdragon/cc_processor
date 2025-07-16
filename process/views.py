from django.shortcuts import render

from process.models import (
    Project,
    SampleStage,
    Abundance,
    StatisticType,
    Statistic,
)

from process.constants import (
    ABUNDANCES_RAW
)

def index(request):
    # Venn diagram of proteins for each project

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



    # Heatmap of some significant proteins
    accession_numbers = ['Q8TF01', 'Q8TF05', 'Q8TF42', 'Q8TF68', 'Q8TF72']

    sample_stages = SampleStage.objects.filter(project__name="SL").order_by('rank')

    samples = [s.name.removeprefix("RO sample ").replace("sample, ", "") for s in sample_stages]

    readings = []

    for accession_number in accession_numbers:
        abundances = Abundance.objects.filter(
            statistic__statistic_type__name = ABUNDANCES_RAW,
            statistic__protein__accession_number = accession_number,
            statistic__protein__project__name = "SL",
            replicate__mean = True
        ).order_by(
            'sample_stage__rank'
        )

        row = [a.reading for a in abundances]

        readings.append(row)

    return render(request, "index.html", {
        "venn_data": sets,
        "accession_numbers_data": accession_numbers,
        "samples_data": samples,
        "readings_data": readings,
    })




