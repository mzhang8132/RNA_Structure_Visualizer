# RNA_Structure_Visualizer
RNA secondary structure predictor based on the Nussinov algorithm with visualization of the best predicted secondary structure.

# Installation

### Run the program through the terminal
1. Clone the repo
2. Install the neccessary libaries (tqdm, networkx, plotly, matplotlib, and nbformat)
```
pip3 install tqdm
pip3 install networkx
pip3 install plotly
pip3 install matplotlib
pip3 install nbformat 
```

### Run the program interactively through a Jupyter notebook
1. Clone the repo
2. Run each code within the Jupyter notebook sequentially

# Usage (Terminal version)

```
--ref <[0, 7]>      Test/reference sequences from RNA Central and RFam database (0-7)
--seq SEQUENCE      RNA sequence
--struct STRUCTURE  RNA structure in dot-parenthesis notation wrapped in quotes
--limit LIMIT       Limit on the number of optimal solutions to return (smaller numbers run faster and won't crash due to memory usage
--vis               Flag to visualize the best predicted structure
```

# Example

### Utilize reference sequence from RFam database from RNA Central

```
python3 nussinov.py --ref 0
```
```
ID: URS00006B0FB5_357808
Generated 1540 optimal structures.
100%|████████████████████████████████████████████████████████████████████████████████████████| 1540/1540 [00:00<00:00, 329065.57it/s]
Best predicted structure accuracy: 0.7
Average predicted structure accuracy: 0.14894241040344977
Best structure: (.().(((((((..()))))))))
```

### Use your own sequence and reference structure

```
python3 nussinov.py --seq AUGG --struct "().."
```
```
Generated 1 optimal structures.
100%|███████████████████████████████████████████████████████████████████████████████████████████████| 1/1 [00:00<00:00, 29330.80it/s]
Best predicted structure accuracy: 1.0
Average predicted structure accuracy: 1.0
```
