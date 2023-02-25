# README 

- Analysis_jupyter
	- cost_estimate: estimation of the cost of operation based on mzedi
	dump logs
	- coupling: combining of analysis to create final filling rate estimates
	- extraction of bin ID numbers and distances from GIS generated files
	- dump_deliveries_analysis: general analysis of dump logs
	- overlay: comparison between extracted ramps and dump logs
	- ramp_rates: main ramp extraction algorithm

- GIS
	- project_ini: general map used for routing estimates, visualisations etc. 

- MATLAB_optim
	- analysis: script to extract specific numbers when analysing results
	- automate_runs: script for running batch optimizations with different parameters
	and saving the results as a single .mat file
	- cost_analysis: calculating labor and operation costs from results
	- full_size_intra_day.m: main optimization script. Requirements: YALMIP, gurobi,
	several input files, lots of time. 
	- graphs: script for generating different types of plots for analysing results
	- archive: old files
	- figures_base: figures from base scenario used in report
	- results: result files with input and output data organized by date and time

- GIS archive: GIS files with missing points, less information

- testbench
	- used for aggregated analysis of all optimization results