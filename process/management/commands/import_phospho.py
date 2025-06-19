import re
import logging
import json

import pandas as pd
from django.core.management.base import BaseCommand
import utilities_basicReader

from process.models import ColumnName, Project, Protein, Phospho, PhosphoReading
# TODO - delete static_mapping
from static_mapping import time_points, data_files, data_files_datakeys, time_points_mapping

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = "import a spreadsheet"

    def add_arguments(self, parser):
        parser.add_argument(
            "--project",
            required=True,
            help="The name of the project to import",
        )

    def handle(self, *args, **options):
        project_name = options["project"]

        if not project_name:
            raise Exception("Please provide a project name")

        logger.info(f"Importing phospho spreadsheet for {project_name}")

        project = Project.objects.get(name=project_name)

        file_path = f"data/{project.phosphoproteome_file}"

        column_names = ColumnName.objects.filter(replicate__project=project)

        cs = Protein.objects.filter(is_contaminant=True, project=project)
        contaminants = [pr.accession_number for pr in cs]

        prs = Protein.objects.filter(is_contaminant=False, project=project)

        if len(prs) == 0:
            raise Exception("No proteins in DB. Have you run import_proteo?")

        proteins = {pr.accession_number: pr for pr in prs}

        PhosphoReading.objects.filter(phospho__protein__project=project).delete()
        Phospho.objects.filter(protein__project=project).delete()

        cns_by_name = {}

        num = 0

        if project_name == "SL":
            for cn in column_names:
                cn_short = re.search(r'IITI_\d{3}_', cn.name)

                cns_by_name[cn_short.group()] = cn

            file_path = f"data/{project.phosphoproteome_file}"

            excel_file = pd.ExcelFile(file_path)

            df = pd.read_excel(excel_file, sheet_name=0)

            row_no = 0

            for index, row in df.iterrows():
                row_no += 1

                accession_number = row[project.proteome_file_accession_number_column_name]

                print(f"Adding phospho row {row_no} {accession_number}")

                if not proteins.get(accession_number):
                    new_protein = Protein.objects.create(
                        project=project, accession_number=accession_number, is_contaminant=False
                    )
                    proteins[accession_number] = new_protein


                protein, created = Protein.objects.get_or_create(
                    project=project, accession_number=accession_number, is_contaminant=False
                )

                mod = row['Modified.Sequence']

                phospho = Phospho.objects.create(
                    protein = protein, mod = mod
                )

                for col in df.columns:
                    col_short = re.search(r'IITI_\d{3}_', col)

                    if col_short is not None:
                        if cn := cns_by_name.get(col_short.group()):
                            reading = row[col]

                            if reading != reading:
                                reading = None

                            print(f"Adding {phospho} {cn} {reading}")

                            PhosphoReading.objects.create(
                                phospho = phospho, column_name=cn, reading=reading,
                            )


            return
        
        if project_name != "ICR":
            raise Exception(f"Unknown project name {project_name}")

        # TODO - split this into two scripts
        # TODO - at least tidy it up

        cns_by_replicate_and_column_name = {}
        time_course_phospho_full = {}
        time_course_phospho_reps = parseTimeCoursePhosphoProteomics(file_path, contaminants)

        for cn in column_names:
            if not cns_by_replicate_and_column_name.get(cn.replicate.name):
                cns_by_replicate_and_column_name[cn.replicate.name] = {}

            cns_by_replicate_and_column_name[cn.replicate.name][cn.sample_stage.name] = cn

        for uniprot_accession in time_course_phospho_reps:
            if not proteins.get(uniprot_accession):
                new_protein = Protein.objects.create(
                    project=project, accession_number=uniprot_accession, is_contaminant=False
                )
                proteins[uniprot_accession] = new_protein

            logger.info(f"Processing protein: {uniprot_accession}")

            if uniprot_accession not in time_course_phospho_full:
                time_course_phospho_full[uniprot_accession] = {"phosphorylation_abundances": {}}

            for mod_key in time_course_phospho_reps[uniprot_accession]["phosphorylation_abundances"]:
                if (
                    mod_key
                    not in time_course_phospho_full[uniprot_accession][
                        "phosphorylation_abundances"
                    ]
                ):
                    num += 1

                    # TODO - is this really needed?
                    time_course_phospho_full[uniprot_accession][
                        "phosphorylation_abundances"
                    ][mod_key] = True

                    mod = time_course_phospho_reps[uniprot_accession][
                        "phosphorylation_abundances"
                    ][mod_key]
                    phosphosite = mod["phosphorylation site"]

                    raw = {
                        "One": {},
                        "Two": {},
                    }

                    # Combine the abundances for a phosphosite from different modifications
                    for modification in mod["peptide_abundances"]["rep_1"]:
                        abundance_rep_1 = mod["peptide_abundances"]["rep_1"][modification]["abundance_rep_1"]
                        abundance_rep_2 = mod["peptide_abundances"]["rep_2"][modification]["abundance_rep_2"]

                        for k, v in abundance_rep_1.items():
                            cur_ab = 0
                            if k in raw["One"]:
                                cur_ab += raw["One"][k]

                            raw["One"][k] = (cur_ab + abundance_rep_1[k])

                        for k, v in abundance_rep_2.items():
                            cur_ab = 0
                            if k in raw["Two"]:
                                cur_ab += raw["Two"][k]

                            raw["Two"][k] = (cur_ab + abundance_rep_2[k])

                phospho = Phospho.objects.create(
                    protein=proteins[uniprot_accession], mod=mod_key, phosphosite=phosphosite
                )

                for replicate_name in raw.keys():
                    for column_name in raw[replicate_name]:
                        reading = raw[replicate_name][column_name]

                        if reading != reading:
                            reading = None

                        col_obj = cns_by_replicate_and_column_name[replicate_name][column_name]

                        PhosphoReading.objects.create(
                            column_name=col_obj, reading=reading, phospho=phospho
                        )


def findModificationPositionRep1(data_point):
    """
    Creates and returns a dictionary that contain information about the positions in the master protein
    and the modification events.
    Example input : ["P20700 1xPhospho [S393]; Q03252 1xPhospho [S407]", "A0A0B4J2F2 1xPhospho [S575]", "O15085 [1452-1473]"
                                    "Q9NYF8 2xPhospho [S290 S/Y]", "Q9P206 2xPhospho [S971 S979]"]

    Example output: modification_info = {'Q9P206': {'S971': {'event': '2xPhospho', 'aa': 'S', 'Position': '971',
            'likelihood': None, 'phosphosite': 'S971', 'Uniprot_Position': {'start': '-', 'end': '-'}},
            'S979': {'event': '2xPhospho', 'aa': 'S', 'Position': '979', 'likelihood': None, 'phosphosite': 'S979', ...}}}}
    """
    
    modification = data_point["Modifications in Master Proteins"]
    positions = findPositionInMasterProtein(data_point["Positions in Master Proteins"])
    modification_info = {}

    regexes = {
        "accession": r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",
        "ptm": r"\dxPhospho",
        "site": r"\[\S+\s*\S*\]",
    }

    modifications = modification.split("; ")
    accession = ""

    if len(modifications) == 0 or modifications == ['']:
        # modifications = "Q07820 [137-176]"
        modification = data_point["Positions in Master Proteins"]
        modifications = [modification]
    for event in modifications:
        try:
            # Uniprot_accession
            uniprot_accession_pattern = re.compile(regexes["accession"])
            accessions = [
                match.group() for match in uniprot_accession_pattern.finditer(event)
            ]
            if len(accessions) > 0:
                accession = accessions[0]
            if accession not in modification_info:
                modification_info[accession] = {}
            # Site
            site_pattern = re.compile(regexes["site"])
            sites = [match.group() for match in site_pattern.finditer(event)]
            if len(sites) == 0:
                # input Q9Y608 2xPhospho [S328(100); S340(99.2)]
                # After split: ['Q9Y608 2xPhospho [S328(100)', 'S340(99.2)]']
                site_pattern = re.compile(r"\[\S+\s*\S*")
                sites = [match.group() for match in site_pattern.finditer(event)]
                if len(sites) == 0:
                    # 'S340(99.2)]'
                    sites = [event]

            sites = findSite(sites)
            # Add the peptide uniprot position on site name of the no specific phosphosites
            # e.g. T/S --> T/S[123-130]
            for p_site, info in sites.items():
                if p_site.find("/") != -1 or len(p_site) == 1:
                    range_name = (
                        "["
                        + str(positions[accession]["start"])
                        + "-"
                        + str(positions[accession]["end"])
                        + "]"
                    )
                    p_site = p_site + range_name
                if (
                    info["phosphosite"].find("/") != -1
                    or len(info["phosphosite"]) == 1
                ):
                    info["phosphosite"] = info["phosphosite"] + range_name
                    info["position"] = info["position"] + range_name

                if p_site not in modification_info[accession]:
                    modification_info[accession][p_site] = {}
                modification_info[accession][p_site]["phosphosite"] = p_site
                modification_info[accession][p_site]["aa"] = info["aa"]
                modification_info[accession][p_site]["position"] = info["position"]
                modification_info[accession][p_site][
                    "uniprot_position"
                ] = positions[accession]
                modification_info[accession][p_site]["likelihood"] = info[
                    "likelihood"
                ]
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(modification)
            print(exc_type, fname, exc_tb.tb_lineno)
            print(e)

    return modification_info

def findPhosphoAbundance(data_source, data_point, uniprot_accession_key):
    """
    Creates and returns a raw phospho abundance dictionary for one protein.
    Data_points: Directly reading from a file.
    """
    ab = {}
    abundance = {uniprot_accession_key: {}}
    for datakey in data_files_datakeys[data_source]:
        value = data_point[datakey]
        if value == "":
            continue
        ab[datakey] = float(value)

    for time_point in time_points:
        col_name = time_points_mapping[data_source][time_point]
        if col_name in ab:
            abundance[uniprot_accession_key][time_point] = ab[col_name]

    return abundance

def findSite(sites):
    """
    Gets the phosphorylation sites for rep1 (merge file) and returns a dictionary with all the information.
    Possible inputs: P0DJD0 [1207-1239]; P49792 [2198-2230] | P27816 1xPhospho [S280] | A0A0B4J2F2 1xPhospho [S575] |
    O15085 [1452-1473] | Q9NYF8 2xPhospho [S290 S/Y] | Q9P206 2xPhospho [S971 S979] | P24928 [1909-1922]; [1923-1936]
    """
    likelihood = None
    final_sites = {}
    for site in sites:
        site = site.split("; ")
        for p_site in site:
            p_site = p_site.split(" ")
            for pp_site in p_site:
                pp_site = pp_site.strip("[").strip("]")
                temp_site_dic = {}
                temp_site_dic[pp_site] = {}
                # [S575]
                if pp_site.find("-") == -1 and pp_site.find("/") == -1:
                    # print("-----> [Sx]")
                    pp_site_full = pp_site.split("(")
                    pp_site = pp_site_full[0]
                    if len(pp_site_full) > 1:
                        likelihood = pp_site_full[1].strip(")")
                    aa = pp_site[0]
                    position = pp_site[1::]
                    if position == "":
                        position = pp_site
                    uniprot_position = {"start": "-", "end": "-"}
                    final_sites[pp_site] = {
                        "aa": aa,
                        "position": position,
                        "phosphosite": pp_site,
                        "uniprot_position": uniprot_position,
                        "likelihood": likelihood,
                    }
                # [S/Y]
                elif pp_site.find("/") != -1:
                    # print("-----> [S/Y]")
                    aa = pp_site
                    position = pp_site
                    uniprot_position = {"start": "-", "end": "-"}
                    final_sites[pp_site] = {
                        "aa": aa,
                        "position": position,
                        "phosphosite": pp_site,
                        "uniprot_position": uniprot_position,
                        "likelihood": likelihood,
                    }
                # [1452-1473]
                elif pp_site.find("-") != -1:
                    # print("-----> [x-x]")
                    aa = "-"
                    position = pp_site
                    site_range = pp_site.split("-")
                    # uniprot_position = {"start": "-", "end": "-"}
                    uniprot_position = {
                        "start": site_range[0],
                        "end": site_range[1],
                    }
                    final_sites[pp_site] = {
                        "aa": aa,
                        "position": position,
                        "phosphosite": pp_site,
                        "uniprot_position": uniprot_position,
                        "likelihood": likelihood,
                    }

    return final_sites

def findPositionInMasterProtein(data_point):
        """
        Finds the positions of the peptide in the master protein.
        Input: "P16402 [35-47]; P10412 [34-46]; P16403 [34-46]"
        Output: {'P16402': {'start': '35', 'end': '47'}, 'P10412': {'start': '34', 'end': '46'}, 'P16403': {'start': '34', 'end': '46'}}
        """
        position_in_master_protein = data_point
        protein_position_info = {}

        positions = position_in_master_protein.split(";")
        for position in positions:
            position = position.strip(" ")
            position_info = position.split(" ")
            for index, item in enumerate(position_info):
                position_info[index] = item.strip("[").strip("]")
                if len(position_info) > 1:
                    position_ranges = position_info[1].split("-")
                else:
                    position_ranges = position_info[0].split("-")
                protein_position_info[position_info[0]] = {
                    "start": position_ranges[0],
                    "end": position_ranges[1],
                }
        return protein_position_info

def parseTimeCoursePhosphoProteomics(file_path, contaminants):
    phospho_rep = {}

    # TODO - delete utilites_basicReader at some point
    data_points_2 = utilities_basicReader.readTableFile(
        file_path,
        byColumn=False,stripQuotes=True,)

    for key, data_point_2 in data_points_2.items():
        if key == 0:
            continue
        try:
            modification_info = findModificationPositionRep1(data_point_2)
            for uniprot_accession_key in modification_info:
                if (
                    uniprot_accession_key == ""
                    or uniprot_accession_key in contaminants
                ):
                    print("+++++ CONTAMINANT FOUND")
                    print(uniprot_accession_key)
                    continue

                abundance_rep_1 = findPhosphoAbundance(
                    "TimeCourse_Phosphoproteomics_rep_1",
                    data_point_2,uniprot_accession_key,)
        
                abundance_rep_2 = findPhosphoAbundance(
                    "TimeCourse_Phosphoproteomics_rep_2",
                    data_point_2, uniprot_accession_key)

                if (
                    len(abundance_rep_1[uniprot_accession_key]) == 0
                    and len(abundance_rep_2[uniprot_accession_key]) == 0
                ):
                    continue
                
                if uniprot_accession_key not in phospho_rep:
                    phospho_rep[uniprot_accession_key] = {
                        "phosphorylation_abundances": {}
                    }

                for phosphosite in modification_info[uniprot_accession_key]:
                    mod_key = modification_info[uniprot_accession_key][phosphosite][
                        "position"
                    ]
                    modification = data_point_2["Modifications"]
                    phosphosite_info = modification_info[uniprot_accession_key][
                        phosphosite
                    ]

                    if (
                        mod_key
                        not in phospho_rep[uniprot_accession_key][
                            "phosphorylation_abundances"
                        ]
                    ):
                        phospho_rep[uniprot_accession_key][
                            "phosphorylation_abundances"
                        ][mod_key] = {
                            "phosphorylation site": phosphosite_info["phosphosite"],
                            "peptide_abundances": {"rep_1":{}, "rep_2":{}},
                            "position_abundances": {
                                "raw": {},
                                "normalised": {"log2_mean": {}, "min-max": {}, "0-max": {}},
                                "imputed": {},
                            },
                            "confidence": {},
                        }

                    phospho_rep[uniprot_accession_key][
                        "phosphorylation_abundances"
                    ][mod_key]["confidence"]["CR07"] = {
                        "protein_groups": data_point_2["Number of Protein Groups"],
                        "no_proteins": data_point_2["Number of Proteins"],
                        "no_isoforms": data_point_2["Number of Isoforms"],
                        "PSMs": data_point_2["Number of PSMs"],
                        "missed_cleavages": data_point_2[
                            "Number of Missed Cleavages"
                        ],
                        "MHplus": data_point_2["Theo MHplus in Da"],
                        "q_value_HT": data_point_2[
                            "Percolator q-Value by Search Engine Sequest HT"
                        ],
                        "PEP_Sequest": data_point_2[
                            "Percolator PEP by Search Engine Sequest HT"
                        ],
                        "XCorr_HT": data_point_2[
                            "XCorr by Search Engine Sequest HT"
                        ],
                        "confidence_Sequest_HT": data_point_2[
                            "Confidence by Search Engine Sequest HT"
                        ],
                    }
                    
                    if (modification not in phospho_rep[uniprot_accession_key][
                            "phosphorylation_abundances"
                        ][mod_key]["peptide_abundances"]["rep_1"]):

                        phospho_rep[uniprot_accession_key][
                            "phosphorylation_abundances"
                        ][mod_key]["peptide_abundances"]["rep_1"][modification] = {
                            "abundance_rep_1": abundance_rep_1[
                                uniprot_accession_key
                            ],
                            "sequence": data_point_2["Annotated Sequence"],
                            "uniprot_position": phosphosite_info[
                                "uniprot_position"
                            ],
                            "aa": phosphosite_info["aa"],
                            "phosphosite": phosphosite_info["phosphosite"],
                            "likelihood": phosphosite_info["likelihood"],
                            "modification": modification,
                        }

                    if (modification not in phospho_rep[uniprot_accession_key][
                            "phosphorylation_abundances"
                        ][mod_key]["peptide_abundances"]["rep_2"]):
                        phospho_rep[uniprot_accession_key][
                            "phosphorylation_abundances"
                        ][mod_key]["peptide_abundances"]["rep_2"][modification] = {
                            "abundance_rep_2": abundance_rep_2[
                                uniprot_accession_key
                            ],
                            "sequence": data_point_2["Annotated Sequence"],
                            "uniprot_position": phosphosite_info[
                                "uniprot_position"
                            ],
                            "aa": phosphosite_info["aa"],
                            "phosphosite": phosphosite_info["phosphosite"],
                            "likelihood": phosphosite_info["likelihood"],
                            "modification": modification,
                        }

        except Exception as e:
            print("this is a problem")
            print(e)

    return phospho_rep