# -*- coding: utf-8 -*-
"""
Created on Friday, July 8, 2016
@author: brettbarbaro
"""

print "hello"

import csv
#import math
#import numpy
#import json

# import data from csv file - Brett
f = "/Users/Brett/Dev/CSVTORECIPE/glycolysis/glycolysis_missinginfo.csv"
all_data = []
#cell_radius = 0.15  # um
#cell_density = 1.07  # g/cc
#protein_content_fraction = 0.163  # 155.0/950.0#0.163#by weight
#cell_volume = 4.0 * math.pi * (math.pow(cell_radius, 3.0) / 3.0)  # cu um
#cell_mass = cell_volume * cell_density * math.pow(10, -12)  # g
#protein_mass = cell_mass * protein_content_fraction  # g
#total_I = [0, 0, 0, 0] 
#col_id = [0, 1]
#total = [total_I, total_IBAQ, total_LFQ]
#oneMol = [total_I, total_IBAQ, total_LFQ]
with open(f, 'rU') as csvfile:
    spamreader = csv.reader(csvfile)
    for row in spamreader:
        all_data.append(row)

#        if len(all_data) == 1: continue
#        for i in range(3):
#            for j in range(4):
#                intensity = row[col_id[i] + j]
#                if intensity == "NaN":
#                    intensity = 0
#                else:
#                    intensity = intensity.replace(",", ".")
#                    total[i][j] += float(intensity)
#nbMol = []
#avgMols = []
#mWeight = []
#names = []



# fill in missing data
#   Ingredients from proteomics:
#     Needs name (the name in the data you have)
#     Needs unique ID e.g. uniprot or ncbi (synthetic cell or http://wholecelldb.stanford.edu)
#     Needs the copy number of concentration (e.g. from mass spec, or http://wholecelldb.stanford.edu/ or david estimation from his painting...)
#     Fetch the sequence
#     Fetch or calculate  the molecular weight
#     From sequence retrieve a template from the PDB. Keep e-value and score…
#     Retrieve the multimer states for the structure using the PDB (biological assembly).
#     From sequence retrieve the localisation if unknow (e.g. cytoplasm, membrane etc,...)
#       http://www.ncbi.nlm.nih.gov/CBBresearch/Lu/subcellular/
#       http://www.imtech.res.in/raghava/pslpred/
#
#     From sequence and template from the PDB build an homology model (e.g. MODELLER or server license for MODELLER is MODELAJE in upper case)
#     If localisation is inner/outer membrane predict/retrieve the orientation of the structure  using http://opm.phar.umich.edu
#     Once a structure is known,
#       Load and center the atoms to the origin
#       apply k-means clustering based on number of atoms.
#       Generate a coarse molecular surface and save it as collada 1.4
#     Generate a recipe file in json.
#
#     Note:
#     http://wholecellkb.stanford.edu/list/Mgenitalium/ProteinComplex contains the mutlimeric states off all  proteins from mycoplasma genitalium. All entry are defined according the gene name (MG_001->MG_526). However for each gene there is a link to the uniprort.
#     Example code are in the github repo of autopack in the script folder
#     For ncbi access the fasta sequence can be done with
#     http://www.ncbi.nlm.nih.gov/protein/296456015?report=fasta&log$=seqview&format=text#

#   export ingredients as json

for x in range(1, len(all_data)):
    name = all_data[x][0]
    pdb = all_data[x][1]
    if all_data[x][2] == '':
        molarity = '0.0'
    else:
        molarity = all_data[x][2]
    ingredient = open("%s.json" % name,"w")
    
    ingredient.write(str('{\n'))
    ingredient.write(str('    "packingMode": "random",\n'))
    ingredient.write(str('    "partners_position": [],\n'))
    ingredient.write(str('    "weight": 0.20000000000000001,\n'))
    ingredient.write(str('    "color": [\n'))
    ingredient.write(str('        0.9,\n'))
    ingredient.write(str('        1,\n'))
    ingredient.write(str('        0\n'))
    ingredient.write(str('    ],\n'))
    ingredient.write(str('    "results": [],\n'))
    ingredient.write(str('    "radii": [\n'))
    ingredient.write(str('        [\n'))
    ingredient.write(str('            5\n'))
    ingredient.write(str('        ]\n'))
    ingredient.write(str('    ],\n'))
    ingredient.write(str('    "cutoff_boundary": null,\n'))
    ingredient.write(str('    "encapsulatingRadius": 20.052597999572754,\n'))
    ingredient.write(str('    "positions2": null,\n'))
    ingredient.write(str('    "rejectionThreshold": 30,\n'))
    ingredient.write(str('    "isAttractor": false,\n'))
    ingredient.write(str('    "Type": "MultiSphere",\n'))
    ingredient.write(str('    "orientBiasRotRangeMax": -3.1415926535897931,\n'))
    ingredient.write(str('    "gradient": "",\n'))
    ingredient.write(str('    "nbMol": 0,\n'))
    ingredient.write(str('    "jitterMax": [\n'))
    ingredient.write(str('        0.20000000000000001,\n'))
    ingredient.write(str('        0.10000000000000001,\n'))
    ingredient.write(str('        0.20000000000000001\n'))
    ingredient.write(str('    ],\n'))
    ingredient.write(str('    "packingPriority": 0,\n'))
    ingredient.write(str('    "rotAxis": [\n'))
    ingredient.write(str('        0.0,\n'))
    ingredient.write(str('        2.0,\n'))
    ingredient.write(str('        1.0\n'))
    ingredient.write(str('    ],\n'))
    ingredient.write(str('    "overwrite_nbMol_value": 0,\n'))
    ingredient.write(str('    "nbJitter": 5,\n'))
    ingredient.write(str('    "molarity": ' + molarity + ',\n'))
    ingredient.write(str('    "useOrientBias": false,\n'))
    ingredient.write(str('    "rotRange": 6.2831000000000001,\n'))
    ingredient.write(str('    "coordsystem": "left",\n'))
    ingredient.write(str('    "meshFile": "' + name + '.dae",\n'))
    ingredient.write(str('    "orientBiasRotRangeMin": -3.1415926535897931,\n'))
    ingredient.write(str('    "perturbAxisAmplitude": 0.10000000000000001,\n'))
    ingredient.write(str('    "proba_binding": 0.5,\n'))
    ingredient.write(str('    "principalVector": [\n'))
    ingredient.write(str('        0.0,\n'))
    ingredient.write(str('        0.0,\n'))
    ingredient.write(str('        -1.0\n'))
    ingredient.write(str('    ],\n'))
    ingredient.write(str('    "properties": {},\n'))
    ingredient.write(str('    "partners_name": [],\n'))
    ingredient.write(str('    "name": "' + name + '",\n'))
    ingredient.write(str('    "positions": [\n'))
    ingredient.write(str('        [\n'))
    ingredient.write(str('            [\n'))
    ingredient.write(str('                0,\n'))
    ingredient.write(str('                0,\n'))
    ingredient.write(str('                0\n'))
    ingredient.write(str('            ]\n'))
    ingredient.write(str('        ]\n'))
    ingredient.write(str('    ],\n'))
    ingredient.write(str('    "meshName": "' + name + '",\n'))
    ingredient.write(str('    "excluded_partners_name": [],\n'))
    ingredient.write(str('    "placeType": "jitter",\n'))
    ingredient.write(str('    "cutoff_surface": 20.052597999572754,\n'))
    ingredient.write(str('    "proba_not_binding": 0.5,\n'))
    ingredient.write(str('    "source": {\n'))
    ingredient.write(str('        "pdb": null,\n'))
    ingredient.write(str('        "transform": {\n'))
    ingredient.write(str('            "center": true\n'))
    ingredient.write(str('        }\n'))
    ingredient.write(str('    },\n'))
    ingredient.write(str('    "use_mesh_rb": false,\n'))
    ingredient.write(str('    "pdb": "' + pdb + '",\n'))
    ingredient.write(str('    "useRotAxis": 1\n'))
    ingredient.write(str('}'))
    
    ingredient.close()

# export recipe as .json
recipe = open("recipe.json","w")

recipe.write(str('{\n'))
recipe.write(str(' "recipe":{\n'))
recipe.write(str('  "version":"1.0",\n'))
recipe.write(str('  "name":"Glycolysis"\n'))
recipe.write(str(' },\n'))
recipe.write(str(' "options":{\n'))
recipe.write(str('  "cancelDialog":false,\n'))
recipe.write(str('  "_hackFreepts":false,\n'))
recipe.write(str('  "windowsSize":1,\n'))
recipe.write(str('  "use_gradient":false,\n'))
recipe.write(str('  "placeMethod":"jitter",\n'))
recipe.write(str('  "saveResult":false,\n'))
recipe.write(str('  "runTimeDisplay":false,\n'))
recipe.write(str('  "overwritePlaceMethod":true,\n'))
recipe.write(str('  "innerGridMethod":"bhtree",\n'))
recipe.write(str('  "boundingBox":[\n'))
recipe.write(str('   [\n'))
recipe.write(str('    -500,\n'))
recipe.write(str('    -500,\n'))
recipe.write(str('    -50\n'))
recipe.write(str('   ],\n'))
recipe.write(str('   [\n'))
recipe.write(str('    500,\n'))
recipe.write(str('    500,\n'))
recipe.write(str('    50\n'))
recipe.write(str('   ]\n'))
recipe.write(str('  ],\n'))
recipe.write(str('  "gradients":[],\n'))
recipe.write(str('  "smallestProteinSize":25,\n'))
recipe.write(str('  "computeGridParams":true,\n'))
recipe.write(str('  "freePtsUpdateThrehod":0.15,\n'))
recipe.write(str('  "pickWeightedIngr":true,\n'))
recipe.write(str('  "_timer":false,\n'))
recipe.write(str('  "ingrLookForNeighbours":false,\n'))
recipe.write(str('  "pickRandPt":true,\n'))
recipe.write(str('  "largestProteinSize":725,\n'))
recipe.write(str('  "resultfile":"",\n'))
recipe.write(str('  "use_periodicity":false,\n'))
recipe.write(str('  "EnviroOnly":false\n'))
recipe.write(str(' },\n'))
recipe.write(str(' "cytoplasme":{\n'))
recipe.write(str('  "ingredients":{\n'))
    

for x in range(1, len(all_data)):
    if all_data[x][0] != '':
        name = all_data[x][0]
        recipe.write(str('   "' + name + '":{\n'))
        recipe.write(str('    "include":"' + name + '.json",\n'))
        recipe.write(str('    "name":"' + name + '"\n'))
        recipe.write(str('   }'))
        if x < (len(all_data) - 1):
            recipe.write(str(',\n'))

recipe.write(str('\n  }\n'))
recipe.write(str(' }\n'))
recipe.write(str('}'))


recipe.close()



print"done"

