#!/usr/bin/env python3
import os
import json
from tabulate import tabulate

MAIN_DIR = os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/..')
FILE_ASSOCIATIONS_PATH = os.path.join(MAIN_DIR, 'variables/file_associations.json')
DOC_DIR = os.path.join(MAIN_DIR, 'docs')


association_table = []

with open(FILE_ASSOCIATIONS_PATH, 'r') as f:
    associations = json.load(f)

for association in associations:
    name = association['name']
    vars = '\n'.join(association['vars'])
    association_table.append((name, vars))

association_table_str = tabulate(
    association_table,
    tablefmt='grid',
    headers=['File', 'Variables']
)

with open(os.path.join(DOC_DIR, 'associations_table.inc'), 'w') as f:
    f.write(association_table_str)
    f.write('\n')
