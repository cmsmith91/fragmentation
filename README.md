# Spatial effects of fragmentation on biodiversity

This repository contains code for the manuscript *Fragmentation mitigates biodiversity loss immediately after habitat destruction*

## General workflow for analyses

1. Conduct parameter exploration of neutral model using [python_code/neutral_mod-amarel15june2021.py](https://github.com/cmsmith91/fragmentation/blob/main/python_code/neutral_mod-amarel15june2021.py) and [r_scripts/look at coms2.R](https://github.com/cmsmith91/fragmentation/blob/main/r_scripts/look%20at%20coms2.R).
2. Run simulations to generate final communities for analysis using [python_code/neutral_mod-amarel15june2021.py](https://github.com/cmsmith91/fragmentation/blob/main/python_code/neutral_mod-amarel15june2021.py).
3. Find the most and least fragmented landscapes in the northeastern USA using [r_scripts/measure_frag.R](https://github.com/cmsmith91/fragmentation/blob/main/r_scripts/measure_frag.R)
4. Run simulation of forest loss using [r_scripts/simulate forest loss.R](https://github.com/cmsmith91/fragmentation/blob/main/r_scripts/simulate%20forest%20loss.R).
5. Summarize results and make figures for the  manuscript using [results_summary.R](https://github.com/cmsmith91/fragmentation/blob/main/r_scripts/results_summary.R), [r_scripts/plot_edge_area.R](https://github.com/cmsmith91/fragmentation/blob/main/r_scripts/plot_edge_area.R), and [r_scripts/methods figure.R](https://github.com/cmsmith91/fragmentation/blob/main/r_scripts/methods%20figure.R).
