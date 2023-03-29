
	# README

	# Directories & Files

	The repository has the following directory tree.

	.
	│   CITATION.cff
	│   CONDUCT.md
	│   LICENSE.md
	│   README.md
	│
	├───data
	│   │   README.md
	│   │
	│   ├───derived_data
	│   │       dist_matrix_ID_filtered.csv
	│   │       mzedi_entry_skip_data_stats.csv
	│   │       README.md
	│   │       skips_indexed_filling_estimates.csv
	│   │
	│   ├───interm_data
	│   │       dist_matrix_ID.csv
	│   │       dist_matrix_ID_filtered.mat
	│   │       filtered_skip_deliveries_at_Mzedi.csv
	│   │       ramping_data.csv
	│   │       ramping_data_merged.csv
	│   │       ramping_data_no_rm_spike.csv
	│   │       skips_indexing.csv
	│   │
	│   └───raw_data
	│           03_Skip_deliveries_at_Mzedi.csv
	│           public-waste-skips-blantyre-malawi.csv
	│           README.md
	│           Skip-fullness-blantyre-raw.csv
	│
	├───docs
	│   │   Gantt_chart.xlsx
	│   │
	│   ├───references
	│   │       bibtex_blantyre_msw.bib
	│   │       ieee.csl
	│   │       semester_bibtex.bib
	│   │
	│   ├───report
	│   │       final_report.pdf
	│   │       report.docx
	│   │
	│   └───slides
	│           project_slides_to_pdf.pdf
	│
	└───src
		│   README.md
		│   testbench.xlsx
		│
		├───Analysis_jupyter
		│   │   cost_estimate.ipynb
		│   │   coupling.ipynb
		│   │   distance_time_mat.ipynb
		│   │   dump_deliveries_analysis.ipynb
		│   │   overlay.ipynb
		│   │   ramp_rates.ipynb
		│   │
		│   ├───.ipynb_checkpoints
		│   │       algo_1-checkpoint.ipynb
		│   │       cost_estimate-checkpoint.ipynb
		│   │       coupling-checkpoint.ipynb
		│   │       cropped-checkpoint.png
		│   │       distance_time_mat-checkpoint.ipynb
		│   │       dump_deliveries_analysis-checkpoint.ipynb
		│   │       overlay-checkpoint.ipynb
		│   │       Q1_vals-checkpoint.csv
		│   │       ramp_rates-checkpoint.ipynb
		│   │       whole_bangwe-checkpoint.png
		│   │
		│   └───figures
		│       ├───.ipynb_checkpoints
		│       ├───mzedi
		│       │   │   weekly_arrivals.png
		│       │   │
		│       │   └───.ipynb_checkpoints
		│       │           weekly_arrivals-checkpoint.png
		│       │
		│       ├───ramps
		│       │   │   Bangwe_inorganic_1_ramps.png
		│       │   │   Bangwe_inorganic_2_ramps.png
		│       │   │   Bangwe_Organic_1_ramps.png
		│       │   │   Bangwe_Organic_2_ramps.png
		│       │   │   BCA_inorganic_1_ramps.png
		│       │   │   BCA_inorganic_2_ramps.png
		│       │   │   BCA_Organic_1_ramps.png
		│       │   │   BCA_Organic_2_ramps.png
		│       │   │   Chigumula_inorganic_1_ramps.png
		│       │   │   Chigumula_inorganic_2_ramps.png
		│       │   │   Chigumula_Organic_1_ramps.png
		│       │   │   Chigumula_Organic_2_ramps.png
		│       │   │   Naizi_inorganic_1_ramps.png
		│       │   │   Naizi_Organic_1_ramps.png
		│       │
		│       ├───ramp_process_naizi_inorganic_1
		│       │   │   bottoms.png
		│       │   │   ramps.png
		│       │   │   spike.png
		│       │   │   tops.png
		│       │
		│       └───raw
		│           │   Bangwe_inorganic_1_raw.png
		│           │   Bangwe_inorganic_2_raw.png
		│           │   Bangwe_Organic_1_raw.png
		│           │   Bangwe_Organic_2_raw.png
		│           │   BCA_inorganic_1_raw.png
		│           │   BCA_inorganic_2_raw.png
		│           │   BCA_Organic_1_raw.png
		│           │   BCA_Organic_2_raw.png
		│           │   Chigumula_inorganic_1_raw.png
		│           │   Chigumula_inorganic_2_raw.png
		│           │   Chigumula_Organic_1_raw.png
		│           │   Chigumula_Organic_2_raw.png
		│           │   Naizi_inorganic_1_raw.png
		│           │   Naizi_Organic_1_raw.png
		│
		├───Archive
		│   │   algo_1.ipynb
		│   │
		│   └───GIS_broken
		│       └───QGIS_map
		│
		├───GIS
		│   └───QGIS_map
		│       │	
		│       ├───figures
		│       │
		│       └───SRTM_files           
		│
		└───MATLAB_optim
			│   analysis.m
			│   automate_runs.m
			│   cost_analysis.m
			│   full_size_intra_day.m
			│   graphs.m
			│
			├───archive
			│       full_size_archived.m
			│       notes_optim.txt
			│       optim_setup.m
			│
			├───figures_base
			│       base_case_exploitation.png
			│       money.fig
			│       money.png
			│       number_skips_vehicles.fig
			│       number_skips_vehicles.png
			│       schedule.fig
			│       schedule.png
			│
			└───results
				├───22-11-21_output
				│       22-11-21_19-52_output.mat
				│       22-11-21_20-12_output.mat
				│
				├───...

	| name         | description                                                                                                                                                                                                                                                                  |
	|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
	| CITATION.cff | A citation file in the citation file format (.cff). Names of contributors and their ORCID iD are added to this file.                                                                                                                                                         |
	| CONDUCT.md   | A code of conduct for this project.                                                                                                                                                                                                                                          |
	| LICENSE.md   | A license for this project. The CC-BY 4.0 license is chosen as the default. 
	|
	| CITATION.cff | A citation file in the citation file format (.cff). Names of contributors and their ORCID iD are added to this file.                                                                                                                                                         |
	| README.md    | A README.md file compiled from README.qmd with format gfm (GitHub Flavoured Markdown)                                                                                                                                                                                        |
	| README.qmd   | A README.qmd file to write up general information about this project.                                                                                                                                                                                                        |
	| data         | A data directory with sub-directories 
	| docs         | A directory for documents that are generated as part of the project. Two sub-folders provide structure for a report and slides.                                                                                                                                                                                                                                    |
	| src          |  folder for scripts that provide machine-readable instructions for data processing.                                                                                                                                                                                  |




# README

# Directories & Files

The repository has the following directory tree.

.

│   CITATION.cff

│   CONDUCT.md

│   LICENSE.md

│   README.md

│

├───data

│   │   README.md

│   │

│   ├───derived_data

│   │       dist_matrix_ID_filtered.csv

│   │       mzedi_entry_skip_data_stats.csv

│   │       README.md

│   │       skips_indexed_filling_estimates.csv

│   │

│   ├───interm_data

│   │       dist_matrix_ID.csv

│   │       dist_matrix_ID_filtered.mat

│   │       filtered_skip_deliveries_at_Mzedi.csv

│   │       ramping_data.csv

│   │       ramping_data_merged.csv

│   │       ramping_data_no_rm_spike.csv

│   │       skips_indexing.csv

│   │

│   └───raw_data

│           03_Skip_deliveries_at_Mzedi.csv

│           public-waste-skips-blantyre-malawi.csv

│           README.md

│           Skip-fullness-blantyre-raw.csv

│

├───docs

│   │   Gantt_chart.xlsx

│   │

│   ├───references

│   │       bibtex_blantyre_msw.bib

│   │       ieee.csl

│   │       semester_bibtex.bib

│   │

│   ├───report

│   │       final_report.pdf

│   │       report.docx

│   │

│   └───slides

│           project_slides_to_pdf.pdf

│

└───src

    │   README.md

    │   testbench.xlsx

    │

    ├───Analysis_jupyter

    │   │   cost_estimate.ipynb

    │   │   coupling.ipynb

    │   │   distance_time_mat.ipynb

    │   │   dump_deliveries_analysis.ipynb

    │   │   overlay.ipynb

    │   │   ramp_rates.ipynb

    │   │

    │   ├───.ipynb_checkpoints

    │   │       algo_1-checkpoint.ipynb

    │   │       cost_estimate-checkpoint.ipynb

    │   │       coupling-checkpoint.ipynb

    │   │       cropped-checkpoint.png

    │   │       distance_time_mat-checkpoint.ipynb

    │   │       dump_deliveries_analysis-checkpoint.ipynb

    │   │       overlay-checkpoint.ipynb

    │   │       Q1_vals-checkpoint.csv

    │   │       ramp_rates-checkpoint.ipynb

    │   │       whole_bangwe-checkpoint.png

    │   │

    │   └───figures

    │       ├───.ipynb_checkpoints

    │       ├───mzedi

    │       │   │   weekly_arrivals.png

    │       │   │

    │       │   └───.ipynb_checkpoints

    │       │           weekly_arrivals-checkpoint.png

    │       │

    │       ├───ramps

    │       │   │   Bangwe_inorganic_1_ramps.png

    │       │   │   Bangwe_inorganic_2_ramps.png

    │       │   │   Bangwe_Organic_1_ramps.png

    │       │   │   Bangwe_Organic_2_ramps.png

    │       │   │   BCA_inorganic_1_ramps.png

    │       │   │   BCA_inorganic_2_ramps.png

    │       │   │   BCA_Organic_1_ramps.png

    │       │   │   BCA_Organic_2_ramps.png

    │       │   │   Chigumula_inorganic_1_ramps.png

    │       │   │   Chigumula_inorganic_2_ramps.png

    │       │   │   Chigumula_Organic_1_ramps.png

    │       │   │   Chigumula_Organic_2_ramps.png

    │       │   │   Naizi_inorganic_1_ramps.png

    │       │   │   Naizi_Organic_1_ramps.png

    │       │

    │       ├───ramp_process_naizi_inorganic_1

    │       │   │   bottoms.png

    │       │   │   ramps.png

    │       │   │   spike.png

    │       │   │   tops.png

    │       │

    │       └───raw

    │           │   Bangwe_inorganic_1_raw.png

    │           │   Bangwe_inorganic_2_raw.png

    │           │   Bangwe_Organic_1_raw.png

    │           │   Bangwe_Organic_2_raw.png

    │           │   BCA_inorganic_1_raw.png

    │           │   BCA_inorganic_2_raw.png

    │           │   BCA_Organic_1_raw.png

    │           │   BCA_Organic_2_raw.png

    │           │   Chigumula_inorganic_1_raw.png

    │           │   Chigumula_inorganic_2_raw.png

    │           │   Chigumula_Organic_1_raw.png

    │           │   Chigumula_Organic_2_raw.png

    │           │   Naizi_inorganic_1_raw.png

    │           │   Naizi_Organic_1_raw.png

    │

    ├───Archive

    │   │   algo_1.ipynb

    │   │

    │   └───GIS_broken

    │       └───QGIS_map

    │

    ├───GIS

    │   └───QGIS_map

    │       │

    │       ├───figures

    │       │

    │       └───SRTM_files           

    │

    └───MATLAB_optim

        │   analysis.m

        │   automate_runs.m

        │   cost_analysis.m

        │   full_size_intra_day.m

        │   graphs.m

        │

        ├───archive

        │       full_size_archived.m

        │       notes_optim.txt

        │       optim_setup.m

        │

        ├───figures_base

        │       base_case_exploitation.png

        │       money.fig

        │       money.png

        │       number_skips_vehicles.fig

        │       number_skips_vehicles.png

        │       schedule.fig

        │       schedule.png

        │

        └───results

            ├───22-11-21_output

            │       22-11-21_19-52_output.mat

            │       22-11-21_20-12_output.mat

            │

            ├───...

| name         | description                                                                                                                                                                                                                                                                  |

|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|

| CITATION.cff | A citation file in the citation file format (.cff). Names of contributors and their ORCID iD are added to this file.                                                                                                                                                         |

| CONDUCT.md   | A code of conduct for this project.                                                                                                                                                                                                                                          |

| LICENSE.md   | A license for this project. The CC-BY 4.0 license is chosen as the default. 

|

| CITATION.cff | A citation file in the citation file format (.cff). Names of contributors and their ORCID iD are added to this file.                                                                                                                                                         |

| README.md    | A README.md file compiled from README.qmd with format gfm (GitHub Flavoured Markdown)                                                                                                                                                                                        |

| README.qmd   | A README.qmd file to write up general information about this project.                                                                                                                                                                                                        |

| data         | A data directory with sub-directories 

| docs         | A directory for documents that are generated as part of the project. Two sub-folders provide structure for a report and slides.                                                                                                                                                                                                                                    |

| src          |  folder for scripts that provide machine-readable instructions for data processing.                                                                                                                                                                                  |