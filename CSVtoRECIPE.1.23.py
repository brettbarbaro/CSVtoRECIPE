# -*- coding: utf-8 -*-
"""
Created on Friday, July 8, 2016
@author: Brett Barbaro, Jared Truong
updated 20160816
"""

# input: a csv file with first four columnns: handle, pdb, molarity, mw
#         only looks at first four columns - others are ignored
# outputs: ingredient .json files for all of the proteins
# outputs: one .json file for the recipe
# outputs: proxy .dae files (tetrahedra)         
# proper dae files need to be downloaded separately
# when all of this is done, autoPACK can build a model with the recipe .json.
# handles must be in proper format - not sure of criteria
# pdb IDs must be in proper format: four character codes
# pdb ID's work as handles

print "hello"

import csv
# import random  # only for color assignments
import os


#import datetime   # optional for timestamp if desired
#current_time = datetime.datetime.strftime(datetime.datetime.now(), '%H.%M.%S')



# import data from csv file - Brett
cwd = os.getcwd() + os.sep

csvname = "genitalium4_cut21"
f = cwd + csvname + ".csv"

pdbpath = cwd + os.sep + 'pdbs' + os.sep


cluster_radius = 10 # radius of spheres in clustered model (min: 5) - doesn't work if it's too small (e.g. nothing lower than 5 worked when I tried it)
print 'cluster_radius =', str(cluster_radius)


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
with open(f, 'rU') as csvfile: # need to open the file in Universal mode so it can read Mac Excel output .csv
    spamreader = csv.reader(csvfile)
    for row in spamreader:
        all_data.append(row)

headers = {'test' : 'headers test works'}

print headers['test']

for num in range(len(all_data[0])):
    headers[all_data[0][num]] = num

print "fetching PDB files"  # - from fetchPDBFilesb.1

import urllib
import csv

# import data from csv file - Brett

if not os.path.isdir('pdbs'):
    print 'making pdbs directory'
    os.mkdir('pdbs')

pdbpath = cwd + '/pdbs/'

''' this code was written by Jared Truong'''
# given pdb ID and path of folder to store files('C:\\Users\\User\\Desktop\\pdbFiles')
# writes pdb file into given location with file name "pdbid.pdb"
def write_pdbFile(pdbid, pdbpath):
    file = open(str(pdbpath) + str(pdbid) + '.pdb', 'w')
    data = fetch_pdb(pdbid)
    file.write(data)
    file.close()

# given pdb ID
# returns pdb file information
def fetch_pdb(pdbid):
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid.upper()
    return urllib.urlopen(url).read()

# given list of pdbid strings, and path of folder to store files
# writes a pdb file for each given IDs
def savePdbFiles(pdbidList, locationPath):
    for i in range(len(pdbidList)):
        write_pdbFile(pdbidList[i], locationPath)
''' this code was written by Jared Truong'''


all_data = []   # this is the structure that contains all of the data

with open(f, 'rU') as csvfile: # need to open the file in Universal mode so it can read Mac Excel output .csv
    spamreader = csv.reader(csvfile)
    for row in spamreader:
        all_data.append(row)

for x in range(1, len(all_data)):                       # for every protein
    pdbid = all_data[x][headers['PDB']]                 # get pdbid
    if pdbid and not os.path.isfile(pdbpath + pdbid + '.pdb'):   # if it's not already there
        write_pdbFile(pdbid, pdbpath)                   # download it

# fill in missing data
#   Ingredients from proteomics:
#     Needs handle (the handle in the data you have)
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

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 16:32:09 2016

@author: Ludovic Autin, Brett Barbaro, Jared Truong

Headers on .csv files must include "HANDLE", "PDB","MOL", "MW", and "NAME"
Entries without "NAME" or where "MOL" is absent or 0 are skipped.
Entries without "PDB" create spheres with volume based on "MW".

"""

import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBList import PDBList
import os
import math
from scipy.cluster.vq import kmeans

def getRadiusFromMW(mw):
    V = mw * 1.21  # protein density = 1.212 cubic angstroms per dalton (Erickson 2009, Biol Proced Online)
    return math.pow((3*V)/(4*math.pi),1.0/3.0)

def buildProxy(handle, pdbid, mw, cluster_radius, pdbfn, surface=False, overwrite=False): # Ludo's code, modified by Brett
    print 'buildProxy for', handle
    # read the pdb
    # build cluster
    # pass pos/radii
    # check for CA only ?
    overwrite = True
    mesh = []
    #structure_id = pdbid
    center=[0,0,0]
    #name = name.split("-")[0]
#    fname = name.split("-")[0]
    #print 'start fetch'
    #pdbid = fetch.retrieve_pdb_file(structure_id) # this method fails about half the time
    #print 'end fetch'
    
#    if len(pdbid) == 4 :
#        if surface :
#            pdbid = getOPM(pdbid)       
#            if pdbid == "":
#                temp_pdbid = fetch.retrieve_pdb_file(structure_id)
#                pdbid = computeOPM(temp_pdbid,pdbid)
#        else :
#            pdbid = fetch.retrieve_pdb_file(structure_id)
#        fname = pdbid
#    else :
#        #parse from protein model portal
#        pdbid = queryPMP(pdbid,name)
#        fname = name
#        if surface and pdbid != "":
#           pdbid = computeOPM(pdbid,name)
#        #http://salilab.org/modbase/retrieve/modbase?databaseID=P21812
#    if pdbid == "" :
#        #gather the radius at least
#        pdbid = 'null'
#        R = getRadiusFromMW(mw)
#        print np.array([[0,0,0]])
#        return [[0,0,0]],[[R]],R,[],[],pdbfn, np.array([[0,0,0]])

    if not os.path.isfile(cwd + os.sep + handle + "_cl.indpolvert") or overwrite:
        s = parser.get_structure(pdbid, pdbfn)
        for m in s.get_models():
            atoms_coord = [atom.coord.tolist() for atom in m.get_atoms() if atom.parent.resname != "DUM"] #
            break       # breaks after first entry
        center = np.average(atoms_coord,axis=0)
        atoms_coord_centerd = np.array(atoms_coord) - center
        R = np.linalg.norm(atoms_coord_centerd,axis=1).max() #R is encapsulating radius - just length of longest vector

#        if cluster_radius == 0:
#            centroids = atoms_coord_centerd
#            nProxy = len(atoms_coord_centerd)
#            cluster_radius = 1
#            return centroids.tolist(), [(np.ones(nProxy)*cluster_radius).tolist(),], R, mesh, pdbfn, center, atoms_coord_centerd

        Vproxy = 4*math.pi*(cluster_radius*cluster_radius*cluster_radius)/3.0
        V = len(atoms_coord) * 10.0 * 1.21  # bbhelp - why * 10???
        nProxy = int(round(V/Vproxy))
        print "cluster ", nProxy, len(atoms_coord)
        if nProxy == 0:
            nProxy = 1
        # ~(1.21 x MW) A**3/molecule
        # V=4*math.pi*(r*r*r)/3.0
        # r =  math.pow((3*V)/(4*math.pi),1.0/3.0)
        # r = math.pow(r3,1.0/3.0)
        #cluster
    #	0.73 cm**3/g  x  10**24 A**3/cm**3  x  molecular weight g/mole
    #	--------------------------------------------------------------
    #			 6.02 x 10**23 molecules/mole
        centroids,_ = kmeans(atoms_coord_centerd,nProxy)
        print 'kmeans completed'
        nProxy = len(centroids) # this is so they are equal in number - kmeans sometimes does not produce nProxy centroids, especially when cluster_radius is low (e.g. <5)
        mesh = []#coarseMolSurface(atoms_coord_centerd.tolist(),None)
        #msms ?

        center = center.tolist()
    else :
        centroids = np.loadtxt(cwd+os.sep+handle+"_cl.indpolvert")
        nProxy = len(centroids)
        R = np.linalg.norm(centroids,axis=1).max()+cluster_radius
    return centroids.tolist(), [(np.ones(nProxy)*cluster_radius).tolist(),], R, mesh, pdbfn, center, atoms_coord_centerd
    

def saveDejaVuMesh(handle, faces, vertices): # Ludo's code, modified by Brett NOTE: .indpolvert must have at least three lines, because cellPACK is expecting triangles.
#    if os.path.isfile(cwd + pdbid + '.indpolvert'):
#        return
    np.savetxt(cwd + handle + ".indpolvert", vertices)
    np.savetxt(cwd + handle + ".indpolface", []) # this is a dummy file - autopack needs it to read in DejaVu



pl = PDBList(pdb=pdbpath)
parser = PDBParser(PERMISSIVE=1)

#if not os.path.isdir('dejavus'):
#    print 'making dejavus directory'
#    os.mkdir('dejavus')
#    
#dejavupath = cwd + os.sep + 'dejavus' + os.sep


print 'writing ingredients'

for x in range(1, len(all_data)):  # for every protein
    handle = all_data[x][headers['HANDLE']]  # get handle, pdb, and molarity
    pdb = all_data[x][headers['PDB']]  # and use to make an ingredient .json
    molarity = all_data[x][headers['MOL']]
    mw = float(all_data[x][headers['MW']])
    pdbfn = pdbpath + str(pdb) + '.pdb'
    
    if not pdb:
        R = getRadiusFromMW(mw)
        proxy = [[0,0,0]],[[R]],R,[],[],pdbfn,[[0,0,0],[0,0,0],[0,0,0]]
    elif handle and molarity and float(molarity) > 0: # if there's none of it in the recipe, don't write the ingredient
        proxy = buildProxy(handle, pdb, mw, cluster_radius, pdbfn, surface=False, overwrite=False)
    else:
        continue
    
    if not pdb:
        pdb = 'null'
    else:
        pdb = '\"' + pdb + '\"'

    positions = proxy[0]
    atoms_coord_centerd = proxy[6]
    saveDejaVuMesh(handle, "", atoms_coord_centerd)    
#    saveDejaVuMesh(handle + '_cl', "", positions)      # don't think we need to save DejaVu file for clusters - that info is contained in the ingredient files.

    ingredient = open("%s.json" % handle,"w")
    
    ingredient.write(str('{\n'))
    ingredient.write(str('    "packingMode": "random",\n'))
    ingredient.write(str('    "partners_position": [],\n'))
    ingredient.write(str('    "weight": 0.20000000000000001,\n'))
    ingredient.write(str('    "color": [\n'))
    ingredient.write(str('        0,\n')) # color doesn't matter at this point (20160816) so I'm picking green!
    ingredient.write(str('        1,\n'))
    ingredient.write(str('        0\n'))
    ingredient.write(str('    ],\n'))
    ingredient.write(str('    "results": [],\n'))
    ingredient.write(str('    "radii": ' + str(proxy[1]) + ',\n'))
    ingredient.write(str('    "cutoff_boundary": null,\n'))
    ingredient.write(str('    "encapsulatingRadius": ' + str(proxy[2]) + ',\n'))
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
    ingredient.write(str('    "meshFile": "' + handle + '",\n'))
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
    ingredient.write(str('    "name": "' + handle + '",\n'))
    ingredient.write(str('    "positions": [\n'))
    ingredient.write(str('        ' + str(positions) + '\n'))
    ingredient.write(str('    ],\n'))
    ingredient.write(str('    "meshName": "' + handle + '",\n'))
    ingredient.write(str('    "excluded_partners_name": [],\n'))
    ingredient.write(str('    "placeType": "jitter",\n'))
    ingredient.write(str('    "cutoff_surface": 20.052597999572754,\n'))
    ingredient.write(str('    "proba_not_binding": 0.5,\n'))
    ingredient.write(str('    "source": {\n'))
    ingredient.write(str('        "pdb": ' + pdb + ',\n'))
    ingredient.write(str('        "transform": {\n'))
    ingredient.write(str('            "center": true\n'))
    ingredient.write(str('        }\n'))
    ingredient.write(str('    },\n'))
    ingredient.write(str('    "use_mesh_rb": false,\n'))
    ingredient.write(str('    "pdb": ' + pdb + ',\n'))
    ingredient.write(str('    "useRotAxis": 1\n'))
    ingredient.write(str('}'))
    
    ingredient.close()


# export recipe as .json
print 'writing recipe'

recipe = open("RECIPE-" + csvname + "-CR" + str(cluster_radius) + ".json","w") #create recipe.json file

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
recipe.write(str('    -500\n'))
recipe.write(str('   ],\n'))
recipe.write(str('   [\n'))
recipe.write(str('    500,\n'))
recipe.write(str('    500,\n'))
recipe.write(str('    500\n'))
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
    

for x in range(1, len(all_data)):           #for each entry in the input .csv, make an entry in the recipe file
    mw = all_data[x][headers['MW']]
    handle = all_data[x][headers['HANDLE']]
    molarity = all_data[x][headers['MOL']]
    if handle and molarity and (float(molarity) > 0):        #if there is no handle in the .csv file, or if molarity=0, the ingredient is skipped - if there is no molarity listed, then an ingredient is created whose molarity can be adjusted later
        if x != 1:         # the first won't be preceeded by a comma. (the last can't have a comma)
            recipe.write(str(',\n'))
        recipe.write(str('   "' + handle + '":{\n'))
        recipe.write(str('    "include":"' + handle + '.json",\n'))
        recipe.write(str('    "name":"' + handle + '"\n'))
        recipe.write(str('   }'))

recipe.write(str('\n  }\n'))
recipe.write(str(' }\n'))
recipe.write(str('}'))


recipe.close()

# writes a tetrahedral dae file as a proxy for every handle

#for x in range(1, len(all_data)): # for every protein
#    handle = all_data[x][headers['HANDLE']]         # get handle and make .dae file
#
#    if all_data[x][0] != '':
#        mw = float(all_data[x][3])
#        radius = .0066141 * (mw ** (1. / 3)) # radius of sphere based on mw; density = 1.212 cubic angstroms per dalton (Erickson 2009, Biol Proced Online)
#
#        dae = open("%s.dae" % handle,"w")
#        
#        dae.write(str('<?xml version="1.0"?>\n'))
#        dae.write(str('<COLLADA xmlns="http://www.collada.org/2005/11/COLLADASchema" version="1.4.1">\n'))
#        dae.write(str('    <asset>\n'))
#        dae.write(str('        <contributor>\n'))
#        dae.write(str('            <authoring_tool>CINEMA4D 17.048 COLLADA Exporter</authoring_tool>\n'))
#        dae.write(str('        </contributor>\n'))
#        dae.write(str('        <created>2016-07-27T01:52:33Z</created>\n'))
#        dae.write(str('        <modified>2016-07-27T01:52:33Z</modified>\n'))
#        dae.write(str('        <unit meter="0.01" name="centimeter"/>\n'))
#        dae.write(str('        <up_axis>Y_UP</up_axis>\n'))
#        dae.write(str('    </asset>\n'))
#        dae.write(str('    <library_geometries>\n'))
#        dae.write(str('        <geometry id="ID4">\n'))
#        dae.write(str('            <mesh>\n'))
#        dae.write(str('                <source id="ID5">\n'))
#        dae.write(str('                    <float_array id="ID6" count="12">' + str(radius * -81.6497) + ' ' + str(radius * -33.3333 ) + ' ' + str(radius * 47.1405 ) + ' ' + str(radius * 81.6497 ) + ' ' + str(radius * -33.3333 ) + ' ' + str(radius * 47.1405 ) + ' ' + str(radius * 0) + ' ' + str(radius * -33.3333) + ' ' + str(radius * -94.2809) + ' ' + str(radius * 0) + ' ' + str(radius * 100) + ' ' + str(radius * -0) + '</float_array>\n'))
#        dae.write(str('                    <technique_common>\n'))
#        dae.write(str('                        <accessor count="4" source="#ID6" stride="3">\n'))
#        dae.write(str('                            <param name="X" type="float"/>\n'))
#        dae.write(str('                            <param name="Y" type="float"/>\n'))
#        dae.write(str('                            <param name="Z" type="float"/>\n'))
#        dae.write(str('                        </accessor>\n'))
#        dae.write(str('                    </technique_common>\n'))
#        dae.write(str('                </source>\n'))
#        dae.write(str('                <source id="ID7">\n'))
#        dae.write(str('                    <float_array id="ID8" count="12">0 -1 -0 0 0.333333 0.942809 0.816497 0.333333 -0.471405 -0.816497 0.333333 -0.471405</float_array>\n'))
#        dae.write(str('                    <technique_common>\n'))
#        dae.write(str('                        <accessor count="4" source="#ID8" stride="3">\n'))
#        dae.write(str('                            <param name="X" type="float"/>\n'))
#        dae.write(str('                            <param name="Y" type="float"/>\n'))
#        dae.write(str('                            <param name="Z" type="float"/>\n'))
#        dae.write(str('                        </accessor>\n'))
#        dae.write(str('                    </technique_common>\n'))
#        dae.write(str('                </source>\n'))
#        dae.write(str('                <source id="ID9">\n'))
#        dae.write(str('                    <float_array id="ID10" count="2">0 1</float_array>\n'))
#        dae.write(str('                    <technique_common>\n'))
#        dae.write(str('                        <accessor count="1" source="#ID10" stride="2">\n'))
#        dae.write(str('                            <param name="S" type="float"/>\n'))
#        dae.write(str('                            <param name="T" type="float"/>\n'))
#        dae.write(str('                        </accessor>\n'))
#        dae.write(str('                    </technique_common>\n'))
#        dae.write(str('                </source>\n'))
#        dae.write(str('                <vertices id="ID11">\n'))
#        dae.write(str('                    <input semantic="POSITION" source="#ID5"/>\n'))
#        dae.write(str('                </vertices>\n'))
#        dae.write(str('                <triangles count="4" material="">\n'))
#        dae.write(str('                    <input offset="0" semantic="VERTEX" source="#ID11"/>\n'))
#        dae.write(str('                    <input offset="1" semantic="NORMAL" source="#ID7"/>\n'))
#        dae.write(str('                    <input offset="2" semantic="TEXCOORD" source="#ID9" set="0"/>\n'))
#        dae.write(str('                    <p>2 0 0 1 0 0 0 0 0 1 1 0 3 1 0 0 1 0 2 2 0 3 2 0 1 2 0 0 3 0 3 3 0 2 3 0</p>\n'))
#        dae.write(str('                </triangles>\n'))
#        dae.write(str('            </mesh>\n'))
#        dae.write(str('        </geometry>\n'))
#        dae.write(str('    </library_geometries>\n'))
#        dae.write(str('    <library_visual_scenes>\n'))
#        dae.write(str('        <visual_scene id="ID1">\n'))
#        dae.write(str('            <node id="ID2" name="' + handle + '">\n'))
#        dae.write(str('                <translate sid="translate">0 0 -0</translate>\n'))
#        dae.write(str('                <rotate sid="rotateY">0 1 0 -0</rotate>\n'))
#        dae.write(str('                <rotate sid="rotateX">1 0 0 0</rotate>\n'))
#        dae.write(str('                <rotate sid="rotateZ">0 0 1 -0</rotate>\n'))
#        dae.write(str('                <scale sid="scale">1 1 1</scale>\n'))
#        dae.write(str('                <node id="ID3" name="' + handle + '">\n'))
#        dae.write(str('                    <translate sid="translate">0 0 -0</translate>\n'))
#        dae.write(str('                    <rotate sid="rotateY">0 1 0 -0</rotate>\n'))
#        dae.write(str('                    <rotate sid="rotateX">1 0 0 0</rotate>\n'))
#        dae.write(str('                    <rotate sid="rotateZ">0 0 1 -0</rotate>\n'))
#        dae.write(str('                    <scale sid="scale">1 1 1</scale>\n'))
#        dae.write(str('                    <instance_geometry url="#ID4"/>\n'))
#        dae.write(str('                </node>\n'))
#        dae.write(str('            </node>\n'))
#        dae.write(str('        </visual_scene>\n'))
#        dae.write(str('    </library_visual_scenes>\n'))
#        dae.write(str('    <scene>\n'))
#        dae.write(str('        <instance_visual_scene url="#ID1"/>\n'))
#        dae.write(str('    </scene>\n'))
#        dae.write(str('</COLLADA>'))
#    
#        dae.close()



print"done"
