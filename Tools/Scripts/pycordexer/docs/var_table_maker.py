#!/usr/bin/env python3
import os
import json
from tabulate import tabulate

MAIN_DIR = os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/..')
CORDEX_VARS_DIR = os.path.join(MAIN_DIR, 'variables/cordex_vars')
DOC_DIR = os.path.join(MAIN_DIR, 'docs')

# Read the variables
CORDEX_VARS = {}
for json_file in os.listdir(CORDEX_VARS_DIR):
    if not json_file.lower().endswith('json'):
        continue

    with open(os.path.join(CORDEX_VARS_DIR, json_file), 'r') as f:
        for var_name, var_description in  json.load(f).items():
            CORDEX_VARS[var_name.lower()] = var_description

var_names = sorted(list(CORDEX_VARS.keys()))
var_descriptions = []

for var_name in var_names:
    var_descriptions.append('')
    actions = CORDEX_VARS[var_name]
    for action in actions:
        for method in action:
            method_name = method[0]
            method_options = method[1]
            if method_name in ('SaveVariableToDisk', 'SaveMultipleVariablesToDisk'):
                if 'new_attributes' in method_options:
                    new_attributes = method_options['new_attributes']
                    if 'long_name' in new_attributes:
                        var_descriptions[-1] = new_attributes['long_name']

var_table = []
for i in range(len(var_descriptions)):
    var_table.append((var_names[i], var_descriptions[i]))

table_str = tabulate(
    var_table,
    tablefmt='grid',
    headers=['Variable', 'Description']
)

with open(os.path.join(DOC_DIR, 'table.inc'), 'w') as f:
    f.write(table_str)
    f.write('\n')
