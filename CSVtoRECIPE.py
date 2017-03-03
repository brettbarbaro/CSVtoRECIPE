# -*- coding: utf-8 -*-
"""
Created on Friday, July 8, 2016
@author: Brett Barbaro, Jared Truong
updated 20160816
"""

# input: a csv file with at least columnns headed: INCLUDE, NAME HANDLE, MOL, PDB
# outputs: ingredient .json files for all of the proteins
# outputs: one .json file for the recipe
# when all of this is done, autoPACK can build a model with the recipe .json.
# handles must be in proper format - not sure of criteria
# PDBIDs must be in proper format: four character codes
# PDBID's work as handles
# "oldnumeric" folder must be copied from
# /Users/mac/Library/Preferences/MAXON/CINEMA 4D R17_89538A46/plugins/ePMV/mgl64/MGLToolsPckgs/numpy" to
# "/Users/mac/anaconda/lib/python2.7/site-packages/numpy" to get the Collada stuff to work

from Bio.PDB.PDBParser import PDBParser
# from Bio.PDB.PDBList import PDBList

import sys

import csv
# import random  # only for color assignments
import os
from collada import Collada
from collada import material
from collada import source
from collada import geometry
from collada import scene
import \
    numpy as np  # "oldnumeric" folder must be copied from /Users/mac/Library/Preferences/MAXON/CINEMA 4D R17_89538A46/plugins/ePMV/mgl64/MGLToolsPckgs/numpy" to "/Users/mac/anaconda/lib/python2.7/site-packages/numpy" to get it to work
import urllib
import math
from scipy.cluster.vq import kmeans

sys.path.insert(0, "/Users/mac/Library/Preferences/MAXON/CINEMA 4D R17_89538A46/plugins/ePMV/mgl64/MGLToolsPckgs/")

print("hello")

# cwd = os.getcwd() + os.sep
model_dir = '/Users/mac/Documents/Alber_model_2/'
csvname = "Alber_model_1A1S_1BXR_1QO1"
recipe_name = model_dir + "RECIPE-" + csvname + ".json"
csvpath = model_dir + csvname + ".csv"

boundingBox = '[[100, 100, 100],[900, 900, 900]]'

surface = False

pdbpath = model_dir + 'PDB' + os.sep

cluster_radius = 'variable'  # radius of spheres in clustered model (min: 5) - doesn't work if it's too small (e.g. nothing lower than 5 worked when I tried it)
print('cluster_radius = ' + str(cluster_radius))


# noinspection PyUnresolvedReferences
def coarseMolSurface(coords, radii, resolution, XYZd=(32, 32, 32), isovalue=6.0, padding=0.0, name='CoarseMolSurface',
                     geom=None):
    print('coarseMolSurface')
    from UTpackages.UTblur import blur
    if radii is None:
        radii = np.ones(len(coords)) * 1.8
    # volarr, origin, span = blur.generateBlurmap(coords, radii, XYZd,resolution, padding = 0.0)
    volarr, origin, span = blur.generateBlurmap(np.ascontiguousarray(coords).tolist(), radii.tolist(), XYZd, resolution,
                                                padding=0.0)
    volarr.shape = (XYZd[0], XYZd[1], XYZd[2])
    volarr = np.ascontiguousarray(np.transpose(volarr), 'f')
    # weights = np.ones(len(radii), "f")
    h = {}
    from Volume.Grid3D import Grid3DF
    maskGrid = Grid3DF(volarr, origin, span, h)
    h['amin'], h['amax'], h['amean'], h['arms'] = maskGrid.stats()
    from UTpackages.UTisocontour import isocontour
    isocontour.setVerboseLevel(0)
    data = maskGrid.data
    origin = np.array(maskGrid.origin).astype('f')
    stepsize = np.array(maskGrid.stepSize).astype('f')
    # add 1 dimension for time steps amd 1 for multiple variables
    if data.dtype.char != np.float32:
        #            print 'converting from ', data.dtype.char
        data = data.astype('f')  # Numeric.Float32)
    newgrid3D = np.ascontiguousarray(np.reshape(np.transpose(data),
                                                (1, 1) + tuple(data.shape)), data.dtype.char)
    #        print "ok"
    ndata = isocontour.newDatasetRegFloat3D(newgrid3D, origin, stepsize)
    #        print "pfff"
    isoc = isocontour.getContour3d(ndata, 0, 0, isovalue,
                                   isocontour.NO_COLOR_VARIABLE)
    vert = np.zeros((isoc.nvert, 3)).astype('f')
    norm = np.zeros((isoc.nvert, 3)).astype('f')
    col = np.zeros(isoc.nvert).astype('f')
    tri = np.zeros((isoc.ntri, 3)).astype('i')
    isocontour.getContour3dData(isoc, vert, norm, col, tri, 0)
    # print vert
    if maskGrid.crystal:
        vert = maskGrid.crystal.toCartesian(vert)
    return vert, norm, tri


def simpleCollada(name, v, f, n, filename):
    print('simpleCollada: ' + name)
    mesh = Collada()
    effect = material.Effect("effect0", [], "phong", diffuse=(1, 0, 0), specular=(0, 1, 0))
    mat = material.Material("material0", "mymaterial", effect)
    mesh.effects.append(effect)
    mesh.materials.append(mat)
    vert_src = source.FloatSource(name + "-v-array", np.array(v).flatten(), ('X', 'Y', 'Z'))
    geom = geometry.Geometry(mesh, "geometry0", name, [vert_src])
    input_list = source.InputList()
    input_list.addInput(0, 'VERTEX', "#" + name + "-v-array")
    triset = geom.createTriangleSet(np.array(f).flatten(), input_list, "materialref")
    geom.primitives.append(triset)
    mesh.geometries.append(geom)
    matnode = scene.MaterialNode("materialref", mat, inputs=[])
    geomnode = scene.GeometryNode(geom, [matnode])
    node = scene.Node(name, children=[geomnode])
    myscene = scene.Scene("myscene", [node])
    mesh.scenes.append(myscene)
    mesh.scene = myscene
    mesh.assetInfo.unitname = "centimeter"
    mesh.assetInfo.unitmeter = 0.01
    mesh.assetInfo.upaxis = "Y_UP"
    mesh.write(filename)


def generateDAEFromPDB(name, coords, filename, resolution):
    print('generateDAEFromPDB: ' + name)
    name = name
    radii = np.array([1.3, ] * len(coords))
    vert, norm, tri = coarseMolSurface(coords, radii, resolution, XYZd=[32, 32, 32], isovalue=6.0, padding=0.0,
                                       name='CoarseMolSurface', geom=None)
    simpleCollada(name, vert, tri, [], filename)


def handleFix(handle):
    handle = handle.replace('.', '_')
    handle = handle.replace(' ', '_')
    handle = handle.replace(',', '_')
    handle = handle.replace('/', '-')
    handle = handle.replace('(', '_')
    handle = handle.replace(')', '_')
    handle = handle.replace('|', '_')
    handle = handle.replace('[', '_')
    handle = handle.replace(']', '_')
    handle = handle.replace(';', '_')
    handle = handle.replace(':',
                            '_')  # have tested characters up to this point - don't work. after here, haven't tested.
    handle = handle.replace('{', '_')
    handle = handle.replace('}', '_')
    handle = handle.replace('\\', '_')
    handle = handle.replace('?', '_')
    handle = handle.replace('<', '_')
    handle = handle.replace('>', '_')
    handle = handle.replace('^', '_')
    handle = handle.replace('%', '_')
    handle = handle.replace('-', '_')  # left these open to add more blacklisted characters later

    x = 1  # this routine adds a number to the handle if it's been used before
    test_handle = handle
    for handles in handle_list:
        if test_handle not in handle_list:
            handle_list.append(test_handle)
            break
        else:
            x += 1
            test_handle = str(handle + '_' + str(x))
    handle = test_handle

    return handle
    # import datetime   # optional for timestamp if desired


# current_time = datetime.datetime.strftime(datetime.datetime.now(), '%H.%M.%S')

# import data from csv file - Brett

all_data = []
handle_list = [None]

# cell_radius = 0.15  # um
# cell_density = 1.07  # g/cc
# protein_content_fraction = 0.163  # 155.0/950.0#0.163#by weight
# cell_volume = 4.0 * math.pi * (math.pow(cell_radius, 3.0) / 3.0)  # cu um
# cell_mass = cell_volume * cell_density * math.pow(10, -12)  # g
# protein_mass = cell_mass * protein_content_fraction  # g

with open(csvpath, 'rU') as csvfile:  # need to open the file in Universal mode so it can read Mac Excel output .csv
    spamreader = csv.reader(csvfile)
    for row in spamreader:
        all_data.append(row)

headers = {'test': 'headers test works'}
print(headers['test'])

for num in range(len(all_data[0])):
    headers[all_data[0][num]] = num  # This establishes a dictionary with the header names in it. After this, columns can be indicated with e.g. "handle = all_data[x][headers['HANDLE']]". The headers must be correctly labeled.

if not os.path.isdir(model_dir + 'PDB'):
    print('making PDB directory')
    os.mkdir(model_dir + 'PDB')


# given PDB ID and path of folder to store files('C:\\Users\\User\\Desktop\\pdbFiles')
# writes PDB file into given location with file name "pdbid.pdb" - this code was written by Jared Truong
def write_pdbFile(pdbid, pdbfilepath):
    print('write_pdbFile')
    data = fetch_pdb(pdbid)
    if data:
        pdbfile = open(str(pdbfilepath) + str(pdbid) + '.pdb', 'w')
        pdbfile.write(data)
        pdbfile.close()
        return True
    else:
        return False


# given PDB ID returns PDB file information - this code was written by Jared Truong
def fetch_pdb(pdbid):
    print('fetching PDB file')
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid.upper()
    return urllib.urlopen(url).read()


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


def getRadiusFromMW(mw):
    mw = float(mw)
    V = mw * 1.212  # protein density = 1.212 cubic angstroms per dalton (Erickson 2009, Biol Proced Online)
    return math.pow((3 * V) / (4 * math.pi), 1.0 / 3.0)


def buildProxy(handle, pdbid, cluster_radius, pdbfn, surface=False, overwrite=False):  # Ludo's code, modified by Brett
    print('buildProxy for ' + handle)
    # read the pdb
    # build cluster
    # pass pos/radii
    # check for CA only ?
    overwrite = True
    mesh = []
    # structure_id = pdbid
    center = [0, 0, 0]
    atoms_coord = []
    atoms_coord_centerd = []
    resolution = -0.1
    xml = None
    # name = name.split("-")[0]
    #    fname = name.split("-")[0]
    # print 'start fetch'
    # pdbid = fetch.retrieve_pdb_file(structure_id) # this method fails about half the time
    # print 'end fetch'

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

    if not os.path.isfile(model_dir + handle + "_cl.indpolvert") or overwrite:
        pdb_struct = parser.get_structure(pdbid, pdbfn)
        for m in pdb_struct.get_models():
            residue_num = 0
            atoms_coord = [atom.coord.tolist() for atom in m.get_atoms() if atom.parent.resname != "DUM"]  #
            atom_num = len(atoms_coord)
            for residues in m.get_residues():
                residue_num += 1
            print('residue_num = ' + str(residue_num))
            print('atoms_coord = ' + str(atom_num))
            if atom_num < 4 * residue_num:
                resolution = -0.05  # changes resolution for alpha-carbon-only structures
            break  # breaks after first entry
        if not atoms_coord:
            print(pdbid, "not found in pdb")
            return None
        center = np.average(atoms_coord, axis=0)  # FIX average position of atoms, weighted by distribution, but not MW
        atoms_coord_centerd = np.array(atoms_coord) - center
        R = np.linalg.norm(atoms_coord_centerd,
                           axis=1).max()  # R is encapsulating radius - just length of longest vector

        #        if cluster_radius == 0:
        #            centroids = atoms_coord_centerd
        #            nProxy = len(atoms_coord_centerd)
        #            cluster_radius = 1
        #            return centroids.tolist(), [(np.ones(nProxy)*cluster_radius).tolist(),], R, mesh, pdbfn, center, atoms_coord_centerd

        V = len(atoms_coord) * 10.0 * 1.21  # bbhelp - why * 10???

        if cluster_radius == 'variable':
            nProxy = 5
            Vproxy = float(V / 5)
            cluster_radius = (((Vproxy * 3) / (math.pi * 4)) ** (
                1. / 3)) * 1.2  # *1.2 is to increase radius a bit so the spheres will touch. there is no support for this.
        else:
            Vproxy = 4 * math.pi * (cluster_radius * cluster_radius * cluster_radius) / 3.0
            nProxy = int(round(V / Vproxy))
            print("clusters:", nProxy, "   atoms:", len(atoms_coord))
        if nProxy == 0:
            nProxy = 1
            # ~(1.21 x MW) A**3/molecule
            # V=4*math.pi*(r*r*r)/3.0
            # r =  math.pow((3*V)/(4*math.pi),1.0/3.0)
            # r = math.pow(r3,1.0/3.0)
            # cluster
            #	0.73 cm**3/g  x  10**24 A**3/cm**3  x  molecular weight g/mole
            #	--------------------------------------------------------------
            #			 6.02 x 10**23 molecules/mole
        centroids, _ = kmeans(atoms_coord_centerd, nProxy)
        print('kmeans completed')
        nProxy = len(
            centroids)  # this is so they are equal in number - kmeans sometimes does not produce nProxy centroids, especially when cluster_radius is low (e.g. <5)
        # mesh = coarseMolSurface(atoms_coord_centerd.tolist(),None)
        # msms ?
        generateDAEFromPDB(handle, atoms_coord_centerd, model_dir + handle + "_coarse.dae", resolution)

        center = center.tolist()
        centroids = centroids.tolist()
    else:
        centroids = np.loadtxt(model_dir + handle + "_cl.indpolvert")
        nProxy = len(centroids)
        R = np.linalg.norm(centroids, axis=1).max() + cluster_radius
        centroids = centroids.tolist()
    return centroids, [(np.ones(nProxy) * cluster_radius).tolist(), ], R, mesh, pdbfn, center, atoms_coord_centerd


def saveDejaVuMesh(handle, faces,
                   vertices):  # Ludo's code, modified by Brett NOTE: .indpolvert must have at least three lines, because cellPACK is expecting triangles.
    #    if os.path.isfile(model_dir + pdbid + '.indpolvert'):
    #        return
    np.savetxt(model_dir + handle + ".indpolvert",
               vertices)  # including all atom positions doesn't significantly slow down results
    np.savetxt(model_dir + handle + ".indpolface", [])  # this is a dummy file - autopack needs it to read in DejaVu


# pl = PDBList(pdb=pdbpath)
parser = PDBParser(PERMISSIVE=True, QUIET=True)  # QUIET=True suppresses warnings about PDB files


# if not os.path.isdir('dejavus'):
#    print 'making dejavus directory'
#    os.mkdir('dejavus')
#
# dejavupath = model_dir + 'dejavus' + os.sep


def writeIngredient(handle, pdb, molarity, mw):
    print('write_Ingredient ' + handle)
    pdbfn = pdbpath + str(pdb) + os.sep + str(pdb) + '.pdb'
    print(pdbfn)
    print(pdb)

    if not pdb:
        R = getRadiusFromMW(mw)
        proxy = [[0, 0, 0]], [[R]], R, [], [], pdbfn, [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        pdb = 'null'
    else:
        proxy = buildProxy(handle, pdb, cluster_radius, pdbfn, surface=False, overwrite=True)

    positions = proxy[0]
    atoms_coord_centerd = proxy[6]
    mesh = proxy[3]
    print(mesh)
    saveDejaVuMesh(handle, "", atoms_coord_centerd)
    #    saveDejaVuMesh(handle + '_cl', "", positions)      # don't think we need to save DejaVu file for clusters - that info is contained in the ingredient files.

    ingredient = open(model_dir + "%s.json" % handle, "w")

    ingredient.write(str('{\n'))
    ingredient.write(str('    "packingMode": "random",\n'))
    ingredient.write(str('    "partners_position": [],\n'))
    ingredient.write(str('    "weight": 0.20000000000000001,\n'))
    ingredient.write(str('    "color": [\n'))
    ingredient.write(str('        0,\n'))  # color doesn't matter at this point (20160816) so I'm picking green!
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
    ingredient.write(str('    "meshFile": "' + handle + "_coarse.dae" + '",\n'))
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
    ingredient.write(str('        "pdb": \"' + pdb + '\",\n'))
    ingredient.write(str('        "transform": {\n'))
    ingredient.write(str('            "center": true\n'))
    ingredient.write(str('        }\n'))
    ingredient.write(str('    },\n'))
    ingredient.write(str('    "use_mesh_rb": false,\n'))
    ingredient.write(str('    "pdb": \"' + pdb + '\",\n'))
    ingredient.write(str('    "useRotAxis": 1\n'))
    ingredient.write(str('}'))

    ingredient.close()


# export recipe as .json
def writeRecipe():
    print('writing recipe')

    recipe = open(recipe_name, "w")  # create recipe.json file

    recipe.write(str('{\n'))
    recipe.write(str(' "recipe":{\n'))
    recipe.write(str('  "version":"1.0",\n'))
    recipe.write(str('  "name":"' + csvname + '"\n'))
    recipe.write(str(' },\n'))
    recipe.write(str(' "options":{\n'))
    recipe.write(str('  "cancelDialog":false,\n'))
    recipe.write(str('  "_hackFreepts":false,\n'))
    recipe.write(str('  "windowsSize":1,\n'))
    recipe.write(str('  "use_gradient":false,\n'))
    recipe.write(str('  "placeMethod":"pandaBullet",\n'))
    recipe.write(str('  "saveResult":false,\n'))
    recipe.write(str('  "runTimeDisplay":false,\n'))
    recipe.write(str('  "overwritePlaceMethod":true,\n'))
    recipe.write(str('  "innerGridMethod":"jordan3",\n'))
    recipe.write(str('  "boundingBox":' + boundingBox + ',\n'))
    recipe.write(str('  "gradients":[],\n'))
    recipe.write(str('  "smallestProteinSize":100,\n'))
    recipe.write(str('  "computeGridParams":true,\n'))
    recipe.write(str('  "freePtsUpdateThrehod":0,\n'))
    recipe.write(str('  "pickWeightedIngr":true,\n'))
    recipe.write(str('  "_timer":false,\n'))
    recipe.write(str('  "ingrLookForNeighbours":false,\n'))
    recipe.write(str('  "pickRandPt":true,\n'))
    recipe.write(str('  "largestProteinSize":100,\n'))
    recipe.write(str('  "resultfile":"",\n'))
    recipe.write(str('  "use_periodicity":false,\n'))
    recipe.write(str('  "EnviroOnly":false\n'))
    recipe.write(str(' },\n'))

    recipe.write(str('''
     "cytoplasme":{
      "ingredients":{
    '''))
    # if surface:
    #     recipe.write(str('''}
    #  },
    #  "compartments":{
    #   "''' + surface + '''":{
    #    "geom":"''' + surface + '''.dae",
    #    "name":"''' + surface + '''",
    #    "surface":{
    #     "ingredients":{}
    #    },
    #    "interior":{
    #     "ingredients":{
    # '''))

    first = True

    for x in range(1, len(all_data)):  # for each entry in the input .csv
        include = all_data[x][headers['INCLUDE']]
        handle = str(all_data[x][headers['HANDLE']])
        pdb = all_data[x][headers['PDB']]
        pdb = pdb.split("_")[0]  # gets rid of any "_X" after PDB name
        molarity = all_data[x][headers['MOL']]
        mw = all_data[x][headers['MW']]
        if include and handle and molarity and (float(molarity) > 0):  # if there is no handle in the .csv file, or if molarity=0, the ingredient is skipped - if there is no molarity listed, then an ingredient is created whose molarity can be adjusted later
            if pdb and not os.path.isdir(model_dir + 'PDB/' + pdb):  # if the PDB's specific directory is not already there
                os.mkdir(model_dir + 'PDB/' + pdb)
                print('making directory for ' + pdb)
                if not os.path.isfile(pdbpath + str(pdb) + os.sep + str(pdb) + '.pdb'):
                    isfile = write_pdbFile(pdb, str(pdbpath + pdb + os.sep))  # download it to its directory
                if not isfile:
                    print(pdb + " NOT FOUND IN PDB")
                    if not mw:
                        print('NO PDB OR MW - SKIPPING ' + pdb)
                    continue
            handle = handleFix(handle)
            if not first:  # the first won't be preceeded by a comma. (the last can't have a comma)
                recipe.write(str(',\n'))
            first = False
            recipe.write(str('     "' + handle + '":{\n'))
            recipe.write(str('      "include":"' + handle + '.json",\n'))
            recipe.write(str('      "name":"' + handle + '"\n'))
            recipe.write(str('     }'))
            writeIngredient(handle, pdb, molarity, mw)

    # if surface:
    #     recipe.write(str('''
    #     }
    #    }'''))

    recipe.write(str('''
      }
     }
    }
    '''))

    recipe.close()


writeRecipe()

os.system('say "Beer time."')

print("done")
