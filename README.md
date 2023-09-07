# Spatial effects of fragmentation on biodiversity

This repository contains code for the manuscript *Fragmentation mitigates biodiversity loss immediately after habitat destruction*

## General workflow for analyses

1. Conduct parameter exploration of neutral model using [python_code/neutral_mod-amarel15june2021.py](https://github.com/cmsmith91/fragmentation/blob/main/python_code/neutral_mod-amarel15june2021.py) and [r_scripts/look at coms2.R](https://github.com/cmsmith91/fragmentation/blob/main/r_scripts/look%20at%20coms2.R).
2. Run simulations to generate final communities for analysis using [python_code/neutral_mod-amarel15june2021.py](https://github.com/cmsmith91/fragmentation/blob/main/python_code/neutral_mod-amarel15june2021.py).
3. Find the most and least fragmented landscapes in the northeastern USA using [r_scripts/measure_frag.R](https://github.com/cmsmith91/fragmentation/blob/main/r_scripts/measure_frag.R)
4. Run simulation of forest loss using [r_scripts/simulate forest loss.R](https://github.com/cmsmith91/fragmentation/blob/main/r_scripts/simulate%20forest%20loss.R).
5. Summarize results and make figures for the  manuscript using [results_summary.R](https://github.com/cmsmith91/fragmentation/blob/main/r_scripts/results_summary.R), [r_scripts/plot_edge_area.R](https://github.com/cmsmith91/fragmentation/blob/main/r_scripts/plot_edge_area.R), and [r_scripts/methods figure.R](https://github.com/cmsmith91/fragmentation/blob/main/r_scripts/methods%20figure.R).


## Description of data in FragData.zip

| File Name  | Description  | 
| :------------ |:---------------| 
| forest_bees.csv      | List of forest specialist bee species in northeastern USA; from [Smith et al 2021](https://www.sciencedirect.com/science/article/pii/S0006320721002548) | 
| forestbee1718_spec.csv      | Data file of all bee specimen records from forest sites in [Smith et al 2021](https://www.sciencedirect.com/science/article/pii/S0006320721002548)        |   
| ForestBeeSiteInfo_final.csv | Data file of site information for forest sites from [Smith et al 2021](https://www.sciencedirect.com/science/article/pii/S0006320721002548)        | 
| hab_maps03dec2021.rds | R list of habitat maps used in the simulations. These were the highest- and lowest-edge landscapes from NLCD data. Each data frame in the list contains the x,y cooordinates of forest and matrix habitat in the landscape        | 
| loss_simulation_n100_output.rds | R list of data frames with output from the forest loss simulations. Each dataframe corresponds to a different landscape size, and each row contains the number of species in the community before and after deforestation   | 
| lu_info03dec2021.rds | R list of data frames with forest landscape data (e.g., percent forest area, edge density) for landscapes we sampled in the norhteastern U.S. Each data frame in the list contains data for a different landscape size | 
| slurm-out | Folder containing text files of raw community data from the neutral model. In the text files, each number is a different species and its position in the vector indicates its x, y coordinate in the 2D map, see [python_code/neutral_mod-amarel15june2021.py](https://github.com/cmsmith91/fragmentation/blob/main/python_code/neutral_mod-amarel15june2021.py)  | 


