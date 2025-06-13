time_points = [
	"Palbo",
	"Late G1_1",
	"G1/S",
	"S",
	"S/G2",
	"G2_2",
	"G2/M_1",
	"M/Early G1",
]

data_files = {
		"TimeCourse_Proteomics_rep_1": "CR07_TMT16plex_FullProt_Proteins.tdt",
		"TimeCourse_Proteomics_rep_2": "CR07_TMT16plex_FullProt_Proteins.tdt",
		"TimeCourse_Phosphoproteomics_rep_1": "20210610_CR07_IMAC_PeptideGroups.tdt",
		"TimeCourse_Phosphoproteomics_rep_2": "20210610_CR07_IMAC_PeptideGroups.tdt",
		"Mitotic_Exit_Proteome":"mitotic_exit_Req93_07JUN2022.tdt",
		"Mitotic_Exit_Phospho":"mitotic_exit_Req93-CR-Davey_IMAC_rep_27JUN22_PeptideGroups.tdt",
		"cached_dataset": "complete_dataset.json",
		"cached_index_protein_names":"total_protein_index.json",
		"cached_TimeCourse_Full_Info": "TimeCourse_Full_info.json",
		"cached_Mitotic_Exit_Full_Info": "Mitotic_Exit_Full_Info.json",
}

data_files_datakeys = {
		"TimeCourse_Proteomics_rep_1": [
			"Abundance: 126, 1, 0",
			"Abundance: 127N, 1, 2",
			"Abundance: 127C, 1, 4",
			"Abundance: 128N, 1, 6",
			"Abundance: 128C, 1, 8",
			"Abundance: 129N, 1, 10",
			"Abundance: 129C, 1, 12",
			"Abundance: 130N, 1, 15"
		],
        
		"TimeCourse_Proteomics_rep_2": [
			"Abundance: 130C, 2, 0",
			"Abundance: 131N, 2, 2",
			"Abundance: 131C, 2, 4",
			"Abundance: 132N, 2, 6",
			"Abundance: 132C, 2, 8",
			"Abundance: 133N, 2, 10",
			"Abundance: 133C, 2, 12",
			"Abundance: 134N, 2, 15"
		],
        
        "TimeCourse_Phosphoproteomics_rep_1": [
			"Abundance F1 126 Sample",
			"Abundance F1 127N Sample",
			"Abundance F1 127C Sample",
			"Abundance F1 128N Sample",
			"Abundance F1 128C Sample",
			"Abundance F1 129N Sample",
			"Abundance F1 129C Sample",
			"Abundance F1 130N Sample"
		],
        
		"TimeCourse_Phosphoproteomics_rep_2": [
			"Abundance F1 130C Sample",
			"Abundance F1 131N Sample",
			"Abundance F1 131C Sample",
			"Abundance F1 132N Sample",
			"Abundance F1 132C Sample",
			"Abundance F1 133N Sample",
			"Abundance F1 133C Sample",
			"Abundance F1 134N Sample"
		],
        
		"Mitotic_Exit_Proteome":[
			"Abundance  127N",
			"Abundance  127C",	
			"Abundance  128N",	
			"Abundance  128C",	
			"Abundance  129N",	
			"Abundance  129C",	
			"Abundance  130N",	
			"Abundance  130C",	
			"Abundance  131N",	
			"Abundance  131C",	
			"Abundance  132N",	
			"Abundance  132C",	
			"Abundance  133N",	
			"Abundance  133C",	
			"Abundance  134N"
		],
	}

time_points_mapping = {
		"TimeCourse_Proteomics_rep_1": {
			"Palbo": "Abundance: 126, 1, 0",
			"Late G1_1": "Abundance: 127N, 1, 2",
			"G1/S": "Abundance: 127C, 1, 4",
			"S": "Abundance: 128N, 1, 6",
			"S/G2": "Abundance: 128C, 1, 8",
			"G2_2": "Abundance: 129N, 1, 10",
			"G2/M_1": "Abundance: 129C, 1, 12",
			"M/Early G1": "Abundance: 130N, 1, 15"
		},
        
		"TimeCourse_Proteomics_rep_2": {
			"Palbo": "Abundance: 130C, 2, 0",
			"Late G1_1": "Abundance: 131N, 2, 2",
			"G1/S": "Abundance: 131C, 2, 4",
			"S": "Abundance: 132N, 2, 6",
			"S/G2": "Abundance: 132C, 2, 8",
			"G2_2": "Abundance: 133N, 2, 10",
			"G2/M_1": "Abundance: 133C, 2, 12",
			"M/Early G1": "Abundance: 134N, 2, 15"
		},
        
		"TimeCourse_Phosphoproteomics_rep_1": {
			"Palbo": "Abundance F1 126 Sample",
			"Late G1_1": "Abundance F1 127N Sample",
			"G1/S": "Abundance F1 127C Sample",
			"S": "Abundance F1 128N Sample",
			"S/G2": "Abundance F1 128C Sample",
			"G2_2": "Abundance F1 129N Sample",
			"G2/M_1": "Abundance F1 129C Sample",
			"M/Early G1": "Abundance F1 130N Sample"
		},
        
		"TimeCourse_Phosphoproteomics_rep_2": {
			"Palbo": "Abundance F1 130C Sample",
			"Late G1_1": "Abundance F1 131N Sample",
			"G1/S": "Abundance F1 131C Sample",
			"S": "Abundance F1 132N Sample",
			"S/G2": "Abundance F1 132C Sample",
			"G2_2": "Abundance F1 133N Sample",
			"G2/M_1": "Abundance F1 133C Sample",
			"M/Early G1": "Abundance F1 134N Sample"
		},

		"Mitotic_Exit_Proteome":{
			"Palbo arrest_R1":"Abundance  127N",
			"Palbo arrest_R2":"Abundance  127C",	
			"Palbo arrest_R3":"Abundance  128N",	
			"DMA arrest_R1":"Abundance  128C",	
			"DMA arrest_R2":"Abundance  129N",	
			"DMA arrest_R3":"Abundance  129C",	
			"DMA release_R1":"Abundance  130N",	
			"DMA release_R2":"Abundance  130C",	
			"DMA release_R3":"Abundance  131N",	
			"Serum starvation arrest_R1":"Abundance  131C",	
			"Serum starvation arrest_R2":"Abundance  132N",	
			"Serum starvation arrest_R3":"Abundance  132C",	
			"Serum starvation release_R1":"Abundance  133N",	
			"Serum starvation release_R2":"Abundance  133C",	
			"Serum starvation release_R3":"Abundance  134N"
		},
	}

replicates_timepoints = {
	"abundance_rep_1" : ["Palbo arrest_R1", "DMA arrest_R1", "DMA release_R1", "Serum starvation arrest_R1", "Serum starvation release_R1"],
	"abundance_rep_2" : ["Palbo arrest_R2", "DMA arrest_R2", "DMA release_R2", "Serum starvation arrest_R2", "Serum starvation release_R2"],
	"abundance_rep_3" : ["Palbo arrest_R3", "DMA arrest_R3", "DMA release_R3", "Serum starvation arrest_R3", "Serum starvation release_R3"]
	}

