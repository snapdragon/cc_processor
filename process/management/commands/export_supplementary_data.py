import logging
from time import sleep
import json
from gprofiler import GProfiler
import pandas as pd
import requests
import io
import csv
from collections import defaultdict

from django.core.management.base import BaseCommand

from process.models import Project, Protein, Phospho, Abundance

from process.constants import ABUNDANCES_RAW

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)

class Command(BaseCommand):
    help = "fetch peptide sequence start positions from uniprot"

    def handle(self, *args, **options):
        logger.info("Export supplementary data for Soliman Labs")

        project = Project.objects.get(name="SL")

        print_percentages(project)

        # _output_CCDs(project)

        # _output_GO_enrichment(project)

        # fetch_mitotic_cell_cycle_proteins(project, "supplementary_data/mitotic_cell_cycle_proteins.csv", True)
        # fetch_mitotic_cell_cycle_proteins(project, "supplementary_data/mitotic_cell_cycle_phospho_proteins.csv", False)

        # _output_well_known_CCDs(project)

        # _output_corum_complexes(project)

        # _output_GO_locations(project)






# List of UniProt accession numbers
accessions = [
    "P06493", "P11802", "P24941", "Q00534", "Q00535", "Q00537", "P14635"
]

# UniProt API endpoint for individual protein entries
def fetch_go_cc_terms(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Error fetching {accession}")
        return []
    data = response.json()
    
    go_terms = []
    for db_ref in data.get("uniProtKBCrossReferences", []):
        if db_ref["database"] == "GO":
            go_id = db_ref["id"]
            term_type = db_ref["properties"][0]["value"]  # F/P/C
            term_name = db_ref["properties"][1]["value"]  # Term name
            if term_type == "C":  # Cellular Component
                go_terms.append((go_id, term_name))
    return go_terms

def _output_corum_complexes(project):
    # TODO - need phospho protein lists too?
    proteins = Protein.objects.filter(
        project=project,
        ccd=True,
    ).exclude(corum_results=None)

    _output_complex(proteins, "supplementary_data/CORUM_CCD_complexes.csv")

    proteins = Protein.objects.filter(
        project=project,
        mitotic_cell_cycle=True,
    ).exclude(corum_results=None)

    _output_complex(proteins, "supplementary_data/CORUM_mitotic_cell_cycle_complexes.csv")


def _output_complex(proteins, filename):
    complexes = {}
    complex_ids = {}

    for protein in proteins:
        for complex in protein.corum_results:
            if complexes.get(complex[0]) is None:
                complexes[complex[0]] = []

            complexes[complex[0]].append(protein.accession_number)
            complex_ids[complex[0]] = complex[1]

    sorted_keys = sorted(complexes, key=lambda k: len(complexes[k]), reverse=True)

    results = []

    for key in sorted_keys:
        results.append({
            "id": key,
            "complex_name": complex_ids[key],
            "count": len(complexes[key]),
            "accession_numbers": complexes[key],
        })

    df = pd.DataFrame(results)

    df.to_csv(filename, index=False)



# Find which proteins and phosphos are associated with GO term 'mitotic cell cycle'
def fetch_mitotic_cell_cycle_proteins(project, filename, is_protein):
    if is_protein:
        proteins = Protein.objects.filter(
            project=project,
            mitotic_cell_cycle=True
        )
    else:
        proteins = Protein.objects.filter(
            project=project,
            phospho_mitotic_cell_cycle=True
        )

    with open(filename, "w") as f:
        f.write("accession_number\n")
        for protein in proteins:
            f.write(f"{protein.accession_number}\n")        

def print_percentages(project):
    # Output percentages of blank cells for proteins and phosphos

    protein_rows = Protein.objects.filter(project=project).count()
    phospho_rows = Phospho.objects.filter(protein__project=project).count()

    raw_protein_abundances_num = Abundance.objects.filter(
        statistic__protein__project = project,
        statistic__statistic_type__name = ABUNDANCES_RAW,
        replicate__mean = False
    ).count()

    raw_phospho_abundances_num = Abundance.objects.filter(
        statistic__phospho__protein__project = project,
        statistic__statistic_type__name = ABUNDANCES_RAW,
        replicate__mean = False
    ).count()

    protein_blank_percentage = round((raw_protein_abundances_num / (protein_rows * 30)) * 100, 2)
    phospho_blank_percentage = round((raw_phospho_abundances_num / (phospho_rows * 30)) * 100, 2)

    print(f"Total populated proteins: {raw_protein_abundances_num} of {protein_rows * 30}: {protein_blank_percentage}% populated")
    print(f"Total populated phosphos: {raw_phospho_abundances_num} of {phospho_rows * 30}: {phospho_blank_percentage}% populated")



def output_CCDs(project):
    # CCD proteins and phosphos

    proteins = Protein.objects.filter(project=project, ccd=True)

    data = []

    # TODO - add more fields to both outputs? See ICR for reference.
    for protein in proteins:
        data.append({"accession_number": protein.accession_number})

    output_file_path = 'supplementary_data/CCDs_protein.csv'

    with open(output_file_path, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=data[0].keys())
        writer.writeheader()
        writer.writerows(data)

    # TODO - should this output phosphopeptided identifiers? If so in
    #   what format?
    # TODO - review all ICR supplementary data and match
    phosphos = Phospho.objects.filter(protein__project=project, ccd=True)

    data = []

    for phospho in phosphos:
        data.append({"accession_number": protein.accession_number})

    output_file_path = 'supplementary_data/CCDs_phospho_protein.csv'

    with open(output_file_path, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=data[0].keys())
        writer.writeheader()
        writer.writerows(data)


def _output_GO_enrichment(project):
    # Output GO enrichment lists

    output_file_path = 'supplementary_data/GO_enrichment_list_protein.csv'

    data = project.protein_go_list

    with open(output_file_path, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=data[0].keys())
        writer.writeheader()
        writer.writerows(data)

    output_file_path = 'supplementary_data/GO_enrichment_list_phospho_protein.csv'

    data = project.phospho_protein_go_list

    with open(output_file_path, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=data[0].keys())
        writer.writeheader()
        writer.writerows(data)

def _output_well_known_CCDs(project):
    # Find which CCDs are in well_known_cell_cycle_regulators.csv
    df = pd.read_csv("data/well_known_cell_cycle_regulators.csv")
    accession_numbers = df['accession_number'].tolist()

    proteins = Protein.objects.filter(
        project=project,
        accession_number__in=accession_numbers,
        ccd=True
    )

    important_accession_numbers = [p.accession_number for p in proteins]

    important_df = df[df['accession_number'].isin(important_accession_numbers)]

    important_df.to_csv('supplementary_data/important_CCDs_found.csv', index=False)

    unimportant_df = df[~df['accession_number'].isin(important_accession_numbers)]

    unimportant_df.to_csv('supplementary_data/important_CCDs_not_found.csv', index=False)


def _output_GO_locations(project):
    # TODO - need phospho protein lists too?
    proteins = Protein.objects.filter(
        project=project,
        ccd=True,
    ).exclude(corum_results=None)

    _output_locations(proteins, "supplementary_data/GO_locations_CCD.csv")

    proteins = Protein.objects.filter(
        project=project,
        mitotic_cell_cycle=True,
    ).exclude(corum_results=None)

    _output_locations(proteins, "supplementary_data/GO_locations_mitotic_cell_cycle.csv")


def _output_locations(proteins, filename):
    locations = {}
    location_ids = {}

    for protein in proteins:
        for location in protein.GO_locations:
            if locations.get(location[0]) is None:
                locations[location[0]] = []

            locations[location[0]].append(protein.accession_number)
            location_ids[location[0]] = location[1]

    sorted_keys = sorted(locations, key=lambda k: len(locations[k]), reverse=True)

    results = []

    for key in sorted_keys:
        results.append({
            "id": key,
            "location_name": location_ids[key],
            "count": len(locations[key]),
            "accession_numbers": locations[key],
        })

    df = pd.DataFrame(results)

    df.to_csv(filename, index=False)
