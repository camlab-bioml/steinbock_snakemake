import argparse
import os
import pandas as pd

# Parser object
parser = argparse.ArgumentParser(prog='Prepare DeepCell', description='Adds deepcell groups for 1 - nuclear and 2 - cytoplasm')

# I/O Arguments
parser.add_argument("-p","--panel", dest = "panel", action = "store", default = "data/data_analysis_training/panel.csv", help = "Name of panel file [default  \"%(default)s\"]")
parser.add_argument("-n","--nuclear", dest = "nuclear", nargs="+", default = ["Ir191"] , help = "Markers for nuclear channels [default  \"%(default)s\"]")
parser.add_argument("-c","--cytoplasm", dest = "cytoplasm", nargs="+", default = ["Er170"], help = "Markers for cytoplasm channels [default  \"%(default)s\"]")
parser.add_argument("-o","--output", dest = "output", action = "store", default = "panel_deepcell.csv", help = "Name of panel output file [default  \"%(default)s\"]")
args = parser.parse_args()

panel = pd.read_csv(args.panel, keep_default_na=False)

# Read in parameters from config
nuclear_channel = args.nuclear
cytoplasm_channel = args.cytoplasm

# Cast deepcell dtype to string
panel.loc[panel['channel'].isin(nuclear_channel), 'deepcell'] = 1
panel.loc[panel['channel'].isin(cytoplasm_channel), 'deepcell'] = 2

panel.to_csv(args.output)

quit()