# -*- coding: utf-8 -*-
"""
Created on Friday, July 8, 2016
@author: Brett Barbaro, Jared Truong
updated 20160816
"""

# input: csv file, top row has headers. Possible headers (must be written exactly as shown):
# INCLUDE - only rows with a value under this header will be processed
# NAME
# MOLARITY
# PDB
# MW
# COLOR
# outputs: ingredient .json files for all of the proteins, one .json file for the recipe, downloads pdb files
# when all of this is done, autoPACK can build a model with the RECIPE_...json file.
# names must be in proper format - not sure of criteria
# PDBIDs must be in proper format: four character codes
# PDBIDs work as names
#
# "oldnumeric" folder must be copied from:
# /Users/mac/Library/Preferences/MAXON/CINEMA 4D R17_89538A46/plugins/ePMV/mgl64/MGLToolsPckgs/numpy/"
# to:
# "/Users/mac/anaconda/lib/python2.7/site-packages/numpy" to get the Collada stuff to work
#
#  The "anaconda" python interpreter package should be installed
#


from Bio.PDB.PDBParser import PDBParser
# from Bio.PDB.PDBList import PDBList

import sys

import csv
import random  # only for color assignments
import os
from collada import Collada
from collada import material
from collada import source
from collada import geometry
from collada import scene
import numpy as np  # "oldnumeric" folder must be copied from /Users/mac/Library/Preferences/MAXON/CINEMA 4D R17_89538A46/plugins/ePMV/mgl64/MGLToolsPckgs/numpy" to "/Users/mac/anaconda/lib/python2.7/site-packages/numpy" to get it to work
import urllib
import math
from scipy.cluster.vq import kmeans

sys.path.insert(0, "/Users/mac/Library/Preferences/MAXON/CINEMA 4D R17_89538A46/plugins/ePMV/mgl64/MGLToolsPckgs/")

print("hello")

csvpath = '/Users/mac/Documents/OLSON/Models/insulin_vesicle/compartment_tests/insulin_vesicle_compartments.csv'
overwrite_ingredients = True
overwrite_dae_files = False

boundingBox = '[[0,0,0],[0,0,0]]'  # When working with dae-defined compartments, this gets adjusted automatically
tree_sphere_radius = 10  # radius of spheres in spheretree clustered model (min: 5) - doesn't work if it's too small (e.g. nothing lower than 5 worked when I tried it); too big is also not good because it becomes low resolution
print('tree_sphere_radius = ' + str(tree_sphere_radius))

# cwd = os.getcwd() + os.sep
head, tail = os.path.split(csvpath)
model_dir = head + '/'
csvname, ext = tail.split('.')
pdbpath = model_dir + 'PDB' + os.sep
recipe_name = model_dir + "RECIPE_" + csvname + ".json"

# current_time = datetime.datetime.strftime(datetime.datetime.now(), '%H.%M.%S')

# import data from csv file - Brett
all_data = []
name_list = [None]

with open(csvpath, 'rU') as csvfile:  # need to open the file in Universal mode so it can read Mac Excel output .csv
    spamreader = csv.reader(csvfile)
    for row in spamreader:
        all_data.append(row)

headers = {'test': 'headers test works'}
print(headers['test'])
for num in range(len(all_data[0])):
    headers[all_data[0][num]] = num  # This establishes a dictionary with the header names in it. After this, columns can be indicated with e.g. "name = all_data[x][headers['NAME']]". The headers must be correctly labeled.

if not os.path.isdir(model_dir + 'PDB'):
    print('making PDB directory')
    os.mkdir(model_dir + 'PDB')


def generateDAEFromPDB(name, coords, filename, resolution):
    print('generateDAEFromPDB: ' + name)
    print('RESOLUTION = ' + str(resolution))
    name = name
    radii = np.array([1.3, ] * len(coords))
    vert, norm, tri = coarseMolSurface(coords, radii, resolution, XYZd=[32, 32, 32], isovalue=6.0, padding=0.0,
                                       name='CoarseMolSurface', geom=None)
    simpleCollada(name, vert, tri, [], filename)


# noinspection PyUnusedLocal
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


# noinspection PyUnresolvedReferences,PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def coarseMolSurface(coords, radii, resolution, XYZd=(32, 32, 32), isovalue=6.0, padding=0.0, name='CoarseMolSurface', geom=None):
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
    newgrid3D = np.ascontiguousarray(np.reshape(np.transpose(data), (1, 1) + tuple(data.shape)), data.dtype.char)
    ndata = isocontour.newDatasetRegFloat3D(newgrid3D, origin, stepsize)
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


# noinspection PyUnusedLocal
def nameFix(name):
    name = name.replace('.', '_')
    name = name.replace(' ', '_')
    name = name.replace(',', '_')
    name = name.replace('/', '-')
    name = name.replace('(', '_')
    name = name.replace(')', '_')
    name = name.replace('|', '_')
    name = name.replace('[', '_')
    name = name.replace(']', '_')
    name = name.replace(';', '_')
    name = name.replace(':', '_')  # characters up to this point don't work. after here, haven't tested.
    name = name.replace('{', '_')
    name = name.replace('}', '_')
    name = name.replace('\\', '_')
    name = name.replace('?', '_')
    name = name.replace('<', '_')
    name = name.replace('>', '_')
    name = name.replace('^', '_')
    name = name.replace('%', '_')
    name = name.replace('-', '_')

    x = 1  # this routine adds a number to the name if it's been used before - e.g. "unknown protein"
    test_name = name
    # noinspection PyUnusedLocal
    for names in name_list:
        if test_name not in name_list:
            name_list.append(test_name)
            break
        else:
            x += 1
            test_name = str(name + '_' + str(x))
    name = test_name

    return name
    # import datetime   # optional for timestamp if desired


# given PDB ID and path of folder to store files('C:\\Users\\User\\Desktop\\pdbFiles'), writes PDB file into given location with file name "pdbid.pdb" - this code was written by Jared Truong
def write_pdbFile(pdbid, pdbfilepath):
    print('write_pdbFile')
    data = fetch_pdb(pdbid)
    print('len(data) = ' + str(len(data)))
    if len(data) > 500:  # 20170801 pdb has started returning short files, ~250 bytes, when a file is not found. This skips those.
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


def getRadiusFromMW(mw):
    mw = float(mw)
    return .66141 * (mw ** (1. / 3))  # radius of minimal sphere based on mw (Erickson 2009, Biol Proced Online)


# noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def buildProxy(name, pdb, cluster_radius, pdbfn, mw, surface=False, overwrite=False):  # Ludo's code, modified by Brett
    print('buildProxy for ' + name)
    overwrite = True
    mesh = []
    center = [0, 0, 0]
    atoms_coord = []
    atoms_coord_centerd = []
    R = 0
    # noinspection PyUnusedLocal
    resolution = -0.1
    # xml = None  # BB not sure what this is for 20170720
    # if not os.path.isfile(model_dir + name + "_cl.indpolvert") or overwrite:
    if not os.path.isfile(pdbfn):
        positions = [[0, 0, 0]]
        R = getRadiusFromMW(mw)
        radii = [[R]]
        if not os.path.isfile(model_dir + name + '_coarse.dae') or overwrite_dae_files is True:
            makeSphericalDae(name, R)
        # proxy = [[0, 0, 0]], [[R]], R, [], [], pdbfn, [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    else:
        pdb_struct = parser.get_structure(pdb, pdbfn)
        atom_multiplier = 1
        for m in pdb_struct.get_models():
            residue_num = 0
            atoms_coord = [atom.coord.tolist() for atom in m.get_atoms() if atom.parent.resname != "DUM"]  #
            if not atoms_coord:
                print(pdb, "not found in pdb")
                return None
            atom_num = len(atoms_coord)
            # noinspection PyUnusedLocal
            for residues in m.get_residues():
                residue_num += 1
            print('residue_num = ' + str(residue_num))
            print('atoms_coord = ' + str(atom_num))
            if atom_num < 4 * residue_num:
                if atom_num == residue_num:
                    resolution = -0.02  # determined by eye
                    atom_multiplier = 7.66  # back-of-envelope calculation using 1b5s.pdb
                    print('ALPHA-CARBON-ONLY DETECTED - RESOLUTION SET TO ' + str(resolution) + ', NUMBER OF ATOMS MULTIPLIED BY ' + str(atom_multiplier))
                else:
                    resolution = -0.06  # THIS IS A GUESS - SHOULD BE FIXED LATER!!!
                    atom_multiplier = 1.55  # back-of-envelope calculation using 1b5s.pdb
                    print('BACKBONE-ONLY DETECTED - RESOLUTION SET TO ' + str(resolution) + ', NUMBER OF ATOMS MULTIPLIED BY ' + str(atom_multiplier))
            break  # breaks after first entry
        center = np.average(atoms_coord, axis=0)  # FIX average position of atoms, weighted by distribution, but not MW
        atoms_coord_centerd = np.array(atoms_coord) - center
        if not os.path.isfile(model_dir + name + '_coarse.dae') or overwrite_dae_files is True:
            generateDAEFromPDB(name, atoms_coord_centerd, model_dir + name + "_coarse.dae", resolution)
        R = np.linalg.norm(atoms_coord_centerd, axis=1).max()  # R is encapsulating radius - just length of longest vector
        V = len(atoms_coord) * 10.0 * 1.21 * atom_multiplier  # bbhelp - why * 10???
        if cluster_radius == 'variable':
            nProxy = 5
            Vproxy = float(V / 5)
            cluster_radius = (((Vproxy * 3) / (math.pi * 4)) ** (1. / 3)) * 1.2  # *1.2 is to increase radius a bit so the spheres will touch. there is no support for this.
        else:
            Vproxy = 4/3.0 * math.pi * (cluster_radius * cluster_radius * cluster_radius)
            nProxy = int(round(V / Vproxy)) + 1  # int() truncates toward 0, +1 prevents nProxy==0, and probably better in general to have more spheres than less
        print("clusters: " + str(nProxy) + "   atoms: " + str(len(atoms_coord)))
        centroids, _ = kmeans(atoms_coord_centerd, nProxy)
        print('kmeans completed')
        nProxy = len(centroids)  # this is so they are equal in number - kmeans sometimes does not produce nProxy centroids, especially when cluster_radius is low (e.g. <5)
        center = center.tolist()
        positions = centroids.tolist()
        radii = [(np.ones(nProxy) * cluster_radius).tolist(), ]

        saveDejaVuMesh(pdb, [], atoms_coord_centerd)  # writes DejaVu mesh containing locations of all atoms

    return positions, radii, R, mesh, pdbfn, center, atoms_coord_centerd


# noinspection PyUnusedLocal
def saveDejaVuMesh(pdb, faces, vertices):  # Ludo's code, modified by Brett NOTE: .indpolvert must have at least three lines, because cellPACK is expecting triangles.
    #    if os.path.isfile(model_dir + pdbid + '.indpolvert'):
    #        return
    np.savetxt(model_dir + 'PDB/' + str(pdb) + os.sep + str(pdb) + ".indpolvert", vertices)  # including all atom positions doesn't significantly slow down results
    np.savetxt(model_dir + 'PDB/' + str(pdb) + os.sep + str(pdb) + ".indpolface", faces)  # this is a dummy file - autopack needs it to read in DejaVu


# writes a tetrahedral dae file as a proxy for every name
def makeTetrahedralDae(name, mw):
    print('makeTetrahedralDae')
    radius = getRadiusFromMW(mw)

    dae = open(model_dir + "%s_coarse.dae" % name, "w")

    dae.write(str('<?xml version="1.0"?>\n'))
    dae.write(str('<COLLADA xmlns="http://www.collada.org/2005/11/COLLADASchema" version="1.4.1">\n'))
    dae.write(str('    <asset>\n'))
    dae.write(str('        <contributor>\n'))
    dae.write(str('            <authoring_tool>CINEMA4D 17.048 COLLADA Exporter</authoring_tool>\n'))
    dae.write(str('        </contributor>\n'))
    dae.write(str('        <created>2016-07-27T01:52:33Z</created>\n'))
    dae.write(str('        <modified>2016-07-27T01:52:33Z</modified>\n'))
    dae.write(str('        <unit meter="0.01" name="centimeter"/>\n'))
    dae.write(str('        <up_axis>Y_UP</up_axis>\n'))
    dae.write(str('    </asset>\n'))
    dae.write(str('    <library_geometries>\n'))
    dae.write(str('        <geometry id="ID4">\n'))
    dae.write(str('            <mesh>\n'))
    dae.write(str('                <source id="ID5">\n'))
    dae.write(str('                    <float_array id="ID6" count="12">' + str(radius * -.816497) + ' ' + str(radius * -.333333) + ' ' + str(radius * .471405) + ' ' + str(radius * .816497) + ' ' + str(radius * -.333333) + ' ' + str(radius * .471405) + ' ' + str(radius * 0) + ' ' + str(radius * -.333333) + ' ' + str(radius * -.942809) + ' ' + str(radius * 0) + ' ' + str(radius * 1) + ' ' + str(radius * 0) + '</float_array>\n'))
    dae.write(str('                    <technique_common>\n'))
    dae.write(str('                        <accessor count="4" source="#ID6" stride="3">\n'))
    dae.write(str('                            <param name="X" type="float"/>\n'))
    dae.write(str('                            <param name="Y" type="float"/>\n'))
    dae.write(str('                            <param name="Z" type="float"/>\n'))
    dae.write(str('                        </accessor>\n'))
    dae.write(str('                    </technique_common>\n'))
    dae.write(str('                </source>\n'))
    dae.write(str('                <source id="ID7">\n'))
    dae.write(str('                    <float_array id="ID8" count="12">0 -1 -0 0 0.333333 0.942809 0.816497 0.333333 -0.471405 -0.816497 0.333333 -0.471405</float_array>\n'))
    dae.write(str('                    <technique_common>\n'))
    dae.write(str('                        <accessor count="4" source="#ID8" stride="3">\n'))
    dae.write(str('                            <param name="X" type="float"/>\n'))
    dae.write(str('                            <param name="Y" type="float"/>\n'))
    dae.write(str('                            <param name="Z" type="float"/>\n'))
    dae.write(str('                        </accessor>\n'))
    dae.write(str('                    </technique_common>\n'))
    dae.write(str('                </source>\n'))
    dae.write(str('                <source id="ID9">\n'))
    dae.write(str('                    <float_array id="ID10" count="2">0 1</float_array>\n'))
    dae.write(str('                    <technique_common>\n'))
    dae.write(str('                        <accessor count="1" source="#ID10" stride="2">\n'))
    dae.write(str('                            <param name="S" type="float"/>\n'))
    dae.write(str('                            <param name="T" type="float"/>\n'))
    dae.write(str('                        </accessor>\n'))
    dae.write(str('                    </technique_common>\n'))
    dae.write(str('                </source>\n'))
    dae.write(str('                <vertices id="ID11">\n'))
    dae.write(str('                    <input semantic="POSITION" source="#ID5"/>\n'))
    dae.write(str('                </vertices>\n'))
    dae.write(str('                <triangles count="4" material="">\n'))
    dae.write(str('                    <input offset="0" semantic="VERTEX" source="#ID11"/>\n'))
    dae.write(str('                    <input offset="1" semantic="NORMAL" source="#ID7"/>\n'))
    dae.write(str('                    <input offset="2" semantic="TEXCOORD" source="#ID9" set="0"/>\n'))
    dae.write(str('                    <p>2 0 0 1 0 0 0 0 0 1 1 0 3 1 0 0 1 0 2 2 0 3 2 0 1 2 0 0 3 0 3 3 0 2 3 0</p>\n'))
    dae.write(str('                </triangles>\n'))
    dae.write(str('            </mesh>\n'))
    dae.write(str('        </geometry>\n'))
    dae.write(str('    </library_geometries>\n'))
    dae.write(str('    <library_visual_scenes>\n'))
    dae.write(str('        <visual_scene id="ID1">\n'))
    dae.write(str('            <node id="ID2" name="' + name + '">\n'))
    dae.write(str('                <translate sid="translate">0 0 -0</translate>\n'))
    dae.write(str('                <rotate sid="rotateY">0 1 0 -0</rotate>\n'))
    dae.write(str('                <rotate sid="rotateX">1 0 0 0</rotate>\n'))
    dae.write(str('                <rotate sid="rotateZ">0 0 1 -0</rotate>\n'))
    dae.write(str('                <scale sid="scale">1 1 1</scale>\n'))
    dae.write(str('                <node id="ID3" name="' + name + '">\n'))
    dae.write(str('                    <translate sid="translate">0 0 -0</translate>\n'))
    dae.write(str('                    <rotate sid="rotateY">0 1 0 -0</rotate>\n'))
    dae.write(str('                    <rotate sid="rotateX">1 0 0 0</rotate>\n'))
    dae.write(str('                    <rotate sid="rotateZ">0 0 1 -0</rotate>\n'))
    dae.write(str('                    <scale sid="scale">1 1 1</scale>\n'))
    dae.write(str('                    <instance_geometry url="#ID4"/>\n'))
    dae.write(str('                </node>\n'))
    dae.write(str('            </node>\n'))
    dae.write(str('        </visual_scene>\n'))
    dae.write(str('    </library_visual_scenes>\n'))
    dae.write(str('    <scene>\n'))
    dae.write(str('        <instance_visual_scene url="#ID1"/>\n'))
    dae.write(str('    </scene>\n'))
    dae.write(str('</COLLADA>'))

    dae.close()


def makeIcosahedralDae(name, mw):
    print('makeIcosahedralDae')
    radius = getRadiusFromMW(mw)

    dae = open(model_dir + "%s_coarse.dae" % name, "w")

    dae.write(str('''<?xml version="1.0"?>
<COLLADA xmlns="http://www.collada.org/2005/11/COLLADASchema" version="1.4.1">
    <asset>
        <contributor>
            <authoring_tool>CINEMA4D 17.055 COLLADA Exporter</authoring_tool>
        </contributor>
        <created>2017-07-28T21:04:06Z</created>
        <modified>2017-07-28T21:04:06Z</modified>
        <unit meter="0.01" name="centimeter"/>
        <up_axis>Y_UP</up_axis>
    </asset>
    <library_geometries>
        <geometry id="ID3">
            <mesh>
                <source id="ID4">
                    <float_array id="ID5" count="36">''' + str(radius * 0) + ' ' + str(radius * 0.346983) + ' ' + str(radius * 0.56143) + ' ' + str(radius * 0) + ' ' + str(radius * -0.346983) + ' ' + str(radius * 0.56143) + ' ' + str(radius * 0.56143) + ' ' + str(radius * 0) + ' ' + str(radius * 0.346983) + ' ' + str(radius * 0.56143) + ' ' + str(radius * 0) + ' ' + str(radius * -0.346983) + ' ' + str(radius * 0) + ' ' + str(radius * 0.346983) + ' ' + str(radius * -0.56143) + ' ' + str(radius * 0) + ' ' + str(radius * -0.346983) + ' ' + str(radius * -0.56143) + ' ' + str(radius * -0.56143) + ' ' + str(radius * 0) + ' ' + str(radius * -0.346983) + ' ' + str(radius * -0.56143) + ' ' + str(radius * 0) + ' ' + str(radius * 0.346983) + ' ' + str(radius * 0.346983) + ' ' + str(radius * 0.56143) + ' ' + str(radius * -0) + ' ' + str(radius * -0.346983) + ' ' + str(radius * 0.56143) + ' ' + str(radius * -0) + ' ' + str(radius * -0.346983) + ' ' + str(radius * -0.56143) + ' ' + str(radius * -0) + ' ' + str(radius * 0.346983) + ' ' + str(radius * -0.56143) + ' ' + str(radius * -0) + '''</float_array>
                    <technique_common>
                        <accessor count="12" source="#ID5" stride="3">
                            <param name="X" type="float"/>
                            <param name="Y" type="float"/>
                            <param name="Z" type="float"/>
                        </accessor>
                    </technique_common>
                </source>
                <source id="ID6">
                    <float_array id="ID7" count="36">-0.525731 0.850651 -0 0.525731 0.850651 -0 0 0.525731 0.850651 0.850651 0 0.525731 0.850651 0 -0.525731 1.50014e-08 0.525731 -0.850651 -0.850651 0 -0.525731 -0.850651 0 0.525731 0.525731 -0.850651 -0 -0.525731 -0.850651 1.50014e-08 0 -0.525731 0.850651 0 -0.525731 -0.850651</float_array>
                    <technique_common>
                        <accessor count="12" source="#ID7" stride="3">
                            <param name="X" type="float"/>
                            <param name="Y" type="float"/>
                            <param name="Z" type="float"/>
                        </accessor>
                    </technique_common>
                </source>
                <source id="ID8">
                    <float_array id="ID9" count="2">0 1</float_array>
                    <technique_common>
                        <accessor count="1" source="#ID9" stride="2">
                            <param name="S" type="float"/>
                            <param name="T" type="float"/>
                        </accessor>
                    </technique_common>
                </source>
                <vertices id="ID10">
                    <input semantic="POSITION" source="#ID4"/>
                </vertices>
                <triangles count="20" material="">
                    <input offset="0" semantic="VERTEX" source="#ID10"/>
                    <input offset="1" semantic="NORMAL" source="#ID6"/>
                    <input offset="2" semantic="TEXCOORD" source="#ID8" set="0"/>
                    <p>0 2 0 8 1 0 9 0 0 0 2 0 2 3 0 8 1 0 2 3 0 3 4 0 8 1 0 3 4 0 4 5 0 8 1 0 4 5 0 9 0 0 8 1 0 6 6 0 9 0 0 4 5 0 7 7 0 9 0 0 6 6 0 7 7 0 0 2 0 9 0 0 1 10 0 10 9 0 11 8 0 1 10 0 11 8 0 2 3 0 11 8 0 3 4 0 2 3 0 11 8 0 5 11 0 3 4 0 11 8 0 10 9 0 5 11 0 10 9 0 6 6 0 5 11 0 10 9 0 7 7 0 6 6 0 10 9 0 1 10 0 7 7 0 0 2 0 7 7 0 1 10 0 0 2 0 1 10 0 2 3 0 3 4 0 5 11 0 4 5 0 5 11 0 6 6 0 4 5 0</p>
                </triangles>
            </mesh>
        </geometry>
    </library_geometries>
    <library_visual_scenes>
        <visual_scene id="ID1">
            <node id="ID2" name="''' + name + '''">
                <translate sid="translate">0 0 -0</translate>
                <rotate sid="rotateY">0 1 0 -0</rotate>
                <rotate sid="rotateX">1 0 0 0</rotate>
                <rotate sid="rotateZ">0 0 1 -0</rotate>
                <scale sid="scale">1 1 1</scale>
                <instance_geometry url="#ID3"/>
            </node>
        </visual_scene>
    </library_visual_scenes>
    <scene>
        <instance_visual_scene url="#ID1"/>
    </scene>
</COLLADA>
'''))

    dae.close()


def makeSphericalDae(name, radius):
    print('makeSphericalDae')

    dae = open(model_dir + '%s_coarse.dae' % name, 'w')

    dae.write('''<?xml version="1.0"?>
<COLLADA xmlns="http://www.collada.org/2005/11/COLLADASchema" version="1.4.1">
    <asset>
        <contributor>
            <authoring_tool>CINEMA4D 17.055 COLLADA Exporter</authoring_tool>
        </contributor>
        <created>2017-07-28T22:34:51Z</created>
        <modified>2017-07-28T22:34:51Z</modified>
        <unit meter="0.01" name="centimeter"/>
        <up_axis>Y_UP</up_axis>
    </asset>
    <library_geometries>
        <geometry id="ID3">
            <mesh>
                <source id="ID4">
                    <float_array id="ID5" count="486">''' + str(radius * 0) + ' ' + str(radius * 0.525731) + ' ' + str(radius * 0.850651) + ' ' + str(radius * 0) + ' ' + str(radius * -0.525731) + ' ' + str(radius * 0.850651) + ' ' + str(radius * 0.850651) + ' ' + str(radius * 0) + ' ' + str(radius * 0.525731) + ' ' + str(radius * 0.850651) + ' ' + str(radius * 0) + ' ' + str(radius * -0.525731) + ' ' + str(radius * 0) + ' ' + str(radius * 0.525731) + ' ' + str(radius * -0.850651) + ' ' + str(radius * 0) + ' ' + str(radius * -0.525731) + ' ' + str(radius * -0.850651) + ' ' + str(radius * -0.850651) + ' ' + str(radius * 0) + ' ' + str(radius * -0.525731) + ' ' + str(radius * -0.850651) + ' ' + str(radius * 0) + ' ' + str(radius * 0.525731) + ' ' + str(radius * 0.525731) + ' ' + str(radius * 0.850651) + ' ' + str(radius * -0) + ' ' + str(radius * -0.525731) + ' ' + str(radius * 0.850651) + ' ' + str(radius * -0) + ' ' + str(radius * -0.525731) + ' ' + str(radius * -0.850651) + ' ' + str(radius * -0) + ' ' + str(radius * 0.525731) + ' ' + str(radius * -0.850651) + ' ' + str(radius * -0) + ' ' + str(radius * 0.295242) + ' ' + str(radius * 0.955423) + ' ' + str(radius * -0) + ' ' + str(radius * 0) + ' ' + str(radius * 1) + ' ' + str(radius * -0) + ' ' + str(radius * -0.295242) + ' ' + str(radius * 0.955423) + ' ' + str(radius * -0) + ' ' + str(radius * 0.147621) + ' ' + str(radius * 0.681718) + ' ' + str(radius * 0.716567) + ' ' + str(radius * 0.309017) + ' ' + str(radius * 0.809017) + ' ' + str(radius * 0.5) + ' ' + str(radius * 0.442863) + ' ' + str(radius * 0.864188) + ' ' + str(radius * 0.238856) + ' ' + str(radius * -0.147621) + ' ' + str(radius * 0.681718) + ' ' + str(radius * 0.716567) + ' ' + str(radius * -0.309017) + ' ' + str(radius * 0.809017) + ' ' + str(radius * 0.5) + ' ' + str(radius * -0.442863) + ' ' + str(radius * 0.864188) + ' ' + str(radius * 0.238856) + ' ' + str(radius * 0.864188) + ' ' + str(radius * 0.238856) + ' ' + str(radius * 0.442863) + ' ' + str(radius * 0.809017) + ' ' + str(radius * 0.5) + ' ' + str(radius * 0.309017) + ' ' + str(radius * 0.681718) + ' ' + str(radius * 0.716567) + ' ' + str(radius * 0.147621) + ' ' + str(radius * 0.238856) + ' ' + str(radius * 0.442863) + ' ' + str(radius * 0.864188) + ' ' + str(radius * 0.5) + ' ' + str(radius * 0.309017) + ' ' + str(radius * 0.809017) + ' ' + str(radius * 0.716567) + ' ' + str(radius * 0.147621) + ' ' + str(radius * 0.681718) + ' ' + str(radius * 0.864188) + ' ' + str(radius * 0.238856) + ' ' + str(radius * -0.442863) + ' ' + str(radius * 0.809017) + ' ' + str(radius * 0.5) + ' ' + str(radius * -0.309017) + ' ' + str(radius * 0.681718) + ' ' + str(radius * 0.716567) + ' ' + str(radius * -0.147621) + ' ' + str(radius * 0.955423) + ' ' + str(radius * 0) + ' ' + str(radius * 0.295242) + ' ' + str(radius * 1) + ' ' + str(radius * 0) + ' ' + str(radius * -0) + ' ' + str(radius * 0.955423) + ' ' + str(radius * 0) + ' ' + str(radius * -0.295242) + ' ' + str(radius * 0.147621) + ' ' + str(radius * 0.681718) + ' ' + str(radius * -0.716567) + ' ' + str(radius * 0.309017) + ' ' + str(radius * 0.809017) + ' ' + str(radius * -0.5) + ' ' + str(radius * 0.442863) + ' ' + str(radius * 0.864188) + ' ' + str(radius * -0.238856) + ' ' + str(radius * 0.716567) + ' ' + str(radius * 0.147621) + ' ' + str(radius * -0.681718) + ' ' + str(radius * 0.5) + ' ' + str(radius * 0.309017) + ' ' + str(radius * -0.809017) + ' ' + str(radius * 0.238856) + ' ' + str(radius * 0.442863) + ' ' + str(radius * -0.864188) + ' ' + str(radius * -0.147621) + ' ' + str(radius * 0.681718) + ' ' + str(radius * -0.716567) + ' ' + str(radius * -0.309017) + ' ' + str(radius * 0.809017) + ' ' + str(radius * -0.5) + ' ' + str(radius * -0.442863) + ' ' + str(radius * 0.864188) + ' ' + str(radius * -0.238856) + ' ' + str(radius * -0.864188) + ' ' + str(radius * 0.238856) + ' ' + str(radius * -0.442863) + ' ' + str(radius * -0.809017) + ' ' + str(radius * 0.5) + ' ' + str(radius * -0.309017) + ' ' + str(radius * -0.681718) + ' ' + str(radius * 0.716567) + ' ' + str(radius * -0.147621) + ' ' + str(radius * -0.238856) + ' ' + str(radius * 0.442863) + ' ' + str(radius * -0.864188) + ' ' + str(radius * -0.5) + ' ' + str(radius * 0.309017) + ' ' + str(radius * -0.809017) + ' ' + str(radius * -0.716567) + ' ' + str(radius * 0.147621) + ' ' + str(radius * -0.681718) + ' ' + str(radius * -0.864188) + ' ' + str(radius * 0.238856) + ' ' + str(radius * 0.442863) + ' ' + str(radius * -0.809017) + ' ' + str(radius * 0.5) + ' ' + str(radius * 0.309017) + ' ' + str(radius * -0.681718) + ' ' + str(radius * 0.716567) + ' ' + str(radius * 0.147621) + ' ' + str(radius * -0.955423) + ' ' + str(radius * 0) + ' ' + str(radius * -0.295242) + ' ' + str(radius * -1) + ' ' + str(radius * 0) + ' ' + str(radius * -0) + ' ' + str(radius * -0.955423) + ' ' + str(radius * 0) + ' ' + str(radius * 0.295242) + ' ' + str(radius * -0.238856) + ' ' + str(radius * 0.442863) + ' ' + str(radius * 0.864188) + ' ' + str(radius * -0.5) + ' ' + str(radius * 0.309017) + ' ' + str(radius * 0.809017) + ' ' + str(radius * -0.716567) + ' ' + str(radius * 0.147621) + ' ' + str(radius * 0.681718) + ' ' + str(radius * -0.295242) + ' ' + str(radius * -0.955423) + ' ' + str(radius * -0) + ' ' + str(radius * 0) + ' ' + str(radius * -1) + ' ' + str(radius * -0) + ' ' + str(radius * 0.295242) + ' ' + str(radius * -0.955423) + ' ' + str(radius * -0) + ' ' + str(radius * -0.147621) + ' ' + str(radius * -0.681718) + ' ' + str(radius * 0.716567) + ' ' + str(radius * -0.309017) + ' ' + str(radius * -0.809017) + ' ' + str(radius * 0.5) + ' ' + str(radius * -0.442863) + ' ' + str(radius * -0.864188) + ' ' + str(radius * 0.238856) + ' ' + str(radius * 0.147621) + ' ' + str(radius * -0.681718) + ' ' + str(radius * 0.716567) + ' ' + str(radius * 0.309017) + ' ' + str(radius * -0.809017) + ' ' + str(radius * 0.5) + ' ' + str(radius * 0.442863) + ' ' + str(radius * -0.864188) + ' ' + str(radius * 0.238856) + ' ' + str(radius * 0.864188) + ' ' + str(radius * -0.238856) + ' ' + str(radius * 0.442863) + ' ' + str(radius * 0.809017) + ' ' + str(radius * -0.5) + ' ' + str(radius * 0.309017) + ' ' + str(radius * 0.681718) + ' ' + str(radius * -0.716567) + ' ' + str(radius * 0.147621) + ' ' + str(radius * 0.238856) + ' ' + str(radius * -0.442863) + ' ' + str(radius * 0.864188) + ' ' + str(radius * 0.5) + ' ' + str(radius * -0.309017) + ' ' + str(radius * 0.809017) + ' ' + str(radius * 0.716567) + ' ' + str(radius * -0.147621) + ' ' + str(radius * 0.681718) + ' ' + str(radius * 0.864188) + ' ' + str(radius * -0.238856) + ' ' + str(radius * -0.442863) + ' ' + str(radius * 0.809017) + ' ' + str(radius * -0.5) + ' ' + str(radius * -0.309017) + ' ' + str(radius * 0.681718) + ' ' + str(radius * -0.716567) + ' ' + str(radius * -0.147621) + ' ' + str(radius * 0.716567) + ' ' + str(radius * -0.147621) + ' ' + str(radius * -0.681718) + ' ' + str(radius * 0.5) + ' ' + str(radius * -0.309017) + ' ' + str(radius * -0.809017) + ' ' + str(radius * 0.238856) + ' ' + str(radius * -0.442863) + ' ' + str(radius * -0.864188) + ' ' + str(radius * 0.147621) + ' ' + str(radius * -0.681718) + ' ' + str(radius * -0.716567) + ' ' + str(radius * 0.309017) + ' ' + str(radius * -0.809017) + ' ' + str(radius * -0.5) + ' ' + str(radius * 0.442863) + ' ' + str(radius * -0.864188) + ' ' + str(radius * -0.238856) + ' ' + str(radius * -0.147621) + ' ' + str(radius * -0.681718) + ' ' + str(radius * -0.716567) + ' ' + str(radius * -0.309017) + ' ' + str(radius * -0.809017) + ' ' + str(radius * -0.5) + ' ' + str(radius * -0.442863) + ' ' + str(radius * -0.864188) + ' ' + str(radius * -0.238856) + ' ' + str(radius * -0.238856) + ' ' + str(radius * -0.442863) + ' ' + str(radius * -0.864188) + ' ' + str(radius * -0.5) + ' ' + str(radius * -0.309017) + ' ' + str(radius * -0.809017) + ' ' + str(radius * -0.716567) + ' ' + str(radius * -0.147621) + ' ' + str(radius * -0.681718) + ' ' + str(radius * -0.864188) + ' ' + str(radius * -0.238856) + ' ' + str(radius * -0.442863) + ' ' + str(radius * -0.809017) + ' ' + str(radius * -0.5) + ' ' + str(radius * -0.309017) + ' ' + str(radius * -0.681718) + ' ' + str(radius * -0.716567) + ' ' + str(radius * -0.147621) + ' ' + str(radius * -0.864188) + ' ' + str(radius * -0.238856) + ' ' + str(radius * 0.442863) + ' ' + str(radius * -0.809017) + ' ' + str(radius * -0.5) + ' ' + str(radius * 0.309017) + ' ' + str(radius * -0.681718) + ' ' + str(radius * -0.716567) + ' ' + str(radius * 0.147621) + ' ' + str(radius * -0.238856) + ' ' + str(radius * -0.442863) + ' ' + str(radius * 0.864188) + ' ' + str(radius * -0.5) + ' ' + str(radius * -0.309017) + ' ' + str(radius * 0.809017) + ' ' + str(radius * -0.716567) + ' ' + str(radius * -0.147621) + ' ' + str(radius * 0.681718) + ' ' + str(radius * 0) + ' ' + str(radius * 0.295242) + ' ' + str(radius * 0.955423) + ' ' + str(radius * 0) + ' ' + str(radius * 0) + ' ' + str(radius * 1) + ' ' + str(radius * 0) + ' ' + str(radius * -0.295242) + ' ' + str(radius * 0.955423) + ' ' + str(radius * 0) + ' ' + str(radius * 0.295242) + ' ' + str(radius * -0.955423) + ' ' + str(radius * 0) + ' ' + str(radius * 0) + ' ' + str(radius * -1) + ' ' + str(radius * 0) + ' ' + str(radius * -0.295242) + ' ' + str(radius * -0.955423) + ' ' + str(radius * -0.16246) + ' ' + str(radius * 0.951057) + ' ' + str(radius * 0.262866) + ' ' + str(radius * 0.16246) + ' ' + str(radius * 0.951057) + ' ' + str(radius * 0.262866) + ' ' + str(radius * 0) + ' ' + str(radius * 0.850651) + ' ' + str(radius * 0.525731) + ' ' + str(radius * 0.587785) + ' ' + str(radius * 0.688191) + ' ' + str(radius * 0.425325) + ' ' + str(radius * 0.688191) + ' ' + str(radius * 0.425325) + ' ' + str(radius * 0.587785) + ' ' + str(radius * 0.425325) + ' ' + str(radius * 0.587785) + ' ' + str(radius * 0.688191) + ' ' + str(radius * 0.850651) + ' ' + str(radius * 0.525731) + ' ' + str(radius * -0) + ' ' + str(radius * 0.951057) + ' ' + str(radius * 0.262866) + ' ' + str(radius * -0.16246) + ' ' + str(radius * 0.951057) + ' ' + str(radius * 0.262866) + ' ' + str(radius * 0.16246) + ' ' + str(radius * 0.587785) + ' ' + str(radius * 0.688191) + ' ' + str(radius * -0.425325) + ' ' + str(radius * 0.425325) + ' ' + str(radius * 0.587785) + ' ' + str(radius * -0.688191) + ' ' + str(radius * 0.688191) + ' ' + str(radius * 0.425325) + ' ' + str(radius * -0.587785) + ' ' + str(radius * 0.16246) + ' ' + str(radius * 0.951057) + ' ' + str(radius * -0.262866) + ' ' + str(radius * -0.16246) + ' ' + str(radius * 0.951057) + ' ' + str(radius * -0.262866) + ' ' + str(radius * 0) + ' ' + str(radius * 0.850651) + ' ' + str(radius * -0.525731) + ' ' + str(radius * -0.425325) + ' ' + str(radius * 0.587785) + ' ' + str(radius * -0.688191) + ' ' + str(radius * -0.587785) + ' ' + str(radius * 0.688191) + ' ' + str(radius * -0.425325) + ' ' + str(radius * -0.688191) + ' ' + str(radius * 0.425325) + ' ' + str(radius * -0.587785) + ' ' + str(radius * -0.951057) + ' ' + str(radius * 0.262866) + ' ' + str(radius * -0.16246) + ' ' + str(radius * -0.850651) + ' ' + str(radius * 0.525731) + ' ' + str(radius * -0) + ' ' + str(radius * -0.951057) + ' ' + str(radius * 0.262866) + ' ' + str(radius * 0.16246) + ' ' + str(radius * -0.587785) + ' ' + str(radius * 0.688191) + ' ' + str(radius * 0.425325) + ' ' + str(radius * -0.425325) + ' ' + str(radius * 0.587785) + ' ' + str(radius * 0.688191) + ' ' + str(radius * -0.688191) + ' ' + str(radius * 0.425325) + ' ' + str(radius * 0.587785) + ' ' + str(radius * 0.16246) + ' ' + str(radius * -0.951057) + ' ' + str(radius * 0.262866) + ' ' + str(radius * -0.16246) + ' ' + str(radius * -0.951057) + ' ' + str(radius * 0.262866) + ' ' + str(radius * 0) + ' ' + str(radius * -0.850651) + ' ' + str(radius * 0.525731) + ' ' + str(radius * 0.688191) + ' ' + str(radius * -0.425325) + ' ' + str(radius * 0.587785) + ' ' + str(radius * 0.587785) + ' ' + str(radius * -0.688191) + ' ' + str(radius * 0.425325) + ' ' + str(radius * 0.425325) + ' ' + str(radius * -0.587785) + ' ' + str(radius * 0.688191) + ' ' + str(radius * 0.951057) + ' ' + str(radius * -0.262866) + ' ' + str(radius * 0.16246) + ' ' + str(radius * 0.951057) + ' ' + str(radius * -0.262866) + ' ' + str(radius * -0.16246) + ' ' + str(radius * 0.850651) + ' ' + str(radius * -0.525731) + ' ' + str(radius * -0) + ' ' + str(radius * 0.688191) + ' ' + str(radius * -0.425325) + ' ' + str(radius * -0.587785) + ' ' + str(radius * 0.425325) + ' ' + str(radius * -0.587785) + ' ' + str(radius * -0.688191) + ' ' + str(radius * 0.587785) + ' ' + str(radius * -0.688191) + ' ' + str(radius * -0.425325) + ' ' + str(radius * 0) + ' ' + str(radius * -0.850651) + ' ' + str(radius * -0.525731) + ' ' + str(radius * -0.16246) + ' ' + str(radius * -0.951057) + ' ' + str(radius * -0.262866) + ' ' + str(radius * 0.16246) + ' ' + str(radius * -0.951057) + ' ' + str(radius * -0.262866) + ' ' + str(radius * -0.425325) + ' ' + str(radius * -0.587785) + ' ' + str(radius * -0.688191) + ' ' + str(radius * -0.688191) + ' ' + str(radius * -0.425325) + ' ' + str(radius * -0.587785) + ' ' + str(radius * -0.587785) + ' ' + str(radius * -0.688191) + ' ' + str(radius * -0.425325) + ' ' + str(radius * -0.951057) + ' ' + str(radius * -0.262866) + ' ' + str(radius * -0.16246) + ' ' + str(radius * -0.951057) + ' ' + str(radius * -0.262866) + ' ' + str(radius * 0.16246) + ' ' + str(radius * -0.850651) + ' ' + str(radius * -0.525731) + ' ' + str(radius * -0) + ' ' + str(radius * -0.688191) + ' ' + str(radius * -0.425325) + ' ' + str(radius * 0.587785) + ' ' + str(radius * -0.425325) + ' ' + str(radius * -0.587785) + ' ' + str(radius * 0.688191) + ' ' + str(radius * -0.587785) + ' ' + str(radius * -0.688191) + ' ' + str(radius * 0.425325) + ' ' + str(radius * -0.262866) + ' ' + str(radius * -0.16246) + ' ' + str(radius * 0.951057) + ' ' + str(radius * -0.525731) + ' ' + str(radius * 0) + ' ' + str(radius * 0.850651) + ' ' + str(radius * -0.262866) + ' ' + str(radius * 0.16246) + ' ' + str(radius * 0.951057) + ' ' + str(radius * 0.525731) + ' ' + str(radius * 0) + ' ' + str(radius * 0.850651) + ' ' + str(radius * 0.262866) + ' ' + str(radius * -0.16246) + ' ' + str(radius * 0.951057) + ' ' + str(radius * 0.262866) + ' ' + str(radius * 0.16246) + ' ' + str(radius * 0.951057) + ' ' + str(radius * 0.262866) + ' ' + str(radius * 0.16246) + ' ' + str(radius * -0.951057) + ' ' + str(radius * 0.262866) + ' ' + str(radius * -0.16246) + ' ' + str(radius * -0.951057) + ' ' + str(radius * 0.525731) + ' ' + str(radius * 0) + ' ' + str(radius * -0.850651) + ' ' + str(radius * -0.262866) + ' ' + str(radius * 0.16246) + ' ' + str(radius * -0.951057) + ' ' + str(radius * -0.525731) + ' ' + str(radius * 0) + ' ' + str(radius * -0.850651) + ' ' + str(radius * -0.262866) + ' ' + str(radius * -0.16246) + ' ' + str(radius * -0.951057) + '''</float_array>
                    <technique_common>
                        <accessor count="162" source="#ID5" stride="3">
                            <param name="X" type="float"/>
                            <param name="Y" type="float"/>
                            <param name="Z" type="float"/>
                        </accessor>
                    </technique_common>
                </source>
                <source id="ID6">
                    <float_array id="ID7" count="486">-0.525731 0.850651 3.01708e-09 -0.292645 0.956221 7.55673e-09 -0.441811 0.864031 0.241356 0 1 -2.52357e-09 -0.154802 0.951018 0.267584 -0.309017 0.809017 0.5 0.292645 0.956221 -5.03782e-09 0.154802 0.951018 0.267584 5.04833e-09 0.855345 0.518058 -0.149166 0.683167 0.714865 0.525731 0.850651 -3.01708e-09 0.441811 0.864031 0.241356 0.309017 0.809017 0.5 0.149166 0.683167 0.714865 3.01708e-09 0.525731 0.850651 0.683167 0.714865 0.149166 0.809017 0.5 0.309017 0.587762 0.683434 0.43296 0.864031 0.241356 0.441811 0.683434 0.43296 0.587762 0.43296 0.587762 0.683434 0.850651 0 0.525731 0.714865 0.149166 0.683167 0.5 0.309017 0.809017 0.241356 0.441811 0.864031 0.683167 0.714865 -0.149166 0.809017 0.5 -0.309017 0.855345 0.518058 2.52416e-09 0.864031 0.241356 -0.441811 0.951018 0.267584 -0.154802 0.951018 0.267584 0.154802 0.850651 3.01708e-09 -0.525731 0.956221 7.55673e-09 -0.292645 1 -2.52357e-09 2.52357e-09 0.956221 0 0.292645 0.441811 0.864031 -0.241356 0.309017 0.809017 -0.5 0.587762 0.683434 -0.43296 0.149166 0.683167 -0.714865 0.43296 0.587762 -0.683434 0.683434 0.43296 -0.587762 -3.01708e-09 0.525731 -0.850651 0.241356 0.441811 -0.864031 0.5 0.309017 -0.809017 0.714865 0.149166 -0.683167 0.154802 0.951018 -0.267584 -0.154802 0.951018 -0.267584 -5.04833e-09 0.855345 -0.518058 -0.441811 0.864031 -0.241356 -0.309017 0.809017 -0.5 -0.149166 0.683167 -0.714865 -0.241356 0.441811 -0.864031 -0.43296 0.587762 -0.683434 -0.5 0.309017 -0.809017 -0.587762 0.683434 -0.43296 -0.683434 0.43296 -0.587762 -0.714865 0.149166 -0.683167 -0.683167 0.714865 -0.149166 -0.809017 0.5 -0.309017 -0.864031 0.241356 -0.441811 -0.850651 0 -0.525731 -0.956221 0 -0.292645 -0.951018 0.267584 -0.154802 -1 -2.52357e-09 -0 -0.855345 0.518058 5.04833e-09 -0.951018 0.267584 0.154802 -0.956221 2.51891e-09 0.292645 -0.683167 0.714865 0.149166 -0.809017 0.5 0.309017 -0.864031 0.241356 0.441811 -0.850651 0 0.525731 -0.587762 0.683434 0.43296 -0.43296 0.587762 0.683434 -0.683434 0.43296 0.587762 -0.241356 0.441811 0.864031 -0.5 0.309017 0.809017 -0.714865 0.149166 0.683167 0.525731 -0.850651 -0 0.292645 -0.956221 5.03782e-09 0.441811 -0.864031 0.241356 0 -1 -2.52357e-09 0.154802 -0.951018 0.267584 0.309017 -0.809017 0.5 -0.292645 -0.956221 -2.51891e-09 -0.154802 -0.951018 0.267584 -5.04833e-09 -0.855345 0.518058 0.149166 -0.683167 0.714865 -0.525731 -0.850651 -6.03417e-09 -0.441811 -0.864031 0.241356 -0.309017 -0.809017 0.5 -0.149166 -0.683167 0.714865 0 -0.525731 0.850651 0.864031 -0.241356 0.441811 0.714865 -0.149166 0.683167 0.809017 -0.5 0.309017 0.683434 -0.43296 0.587762 0.5 -0.309017 0.809017 0.683167 -0.714865 0.149166 0.587762 -0.683434 0.43296 0.43296 -0.587762 0.683434 0.241356 -0.441811 0.864031 0.951018 -0.267584 0.154802 0.951018 -0.267584 -0.154802 0.855345 -0.518058 -5.04833e-09 0.864031 -0.241356 -0.441811 0.809017 -0.5 -0.309017 0.683167 -0.714865 -0.149166 0.714865 -0.149166 -0.683167 0.5 -0.309017 -0.809017 0.683434 -0.43296 -0.587762 0.241356 -0.441811 -0.864031 0.43296 -0.587762 -0.683434 0.587762 -0.683434 -0.43296 0 -0.525731 -0.850651 0.149166 -0.683167 -0.714865 0.309017 -0.809017 -0.5 0.441811 -0.864031 -0.241356 -0.149166 -0.683167 -0.714865 -0.309017 -0.809017 -0.5 2.52416e-09 -0.855345 -0.518058 -0.441811 -0.864031 -0.241356 -0.154802 -0.951018 -0.267584 0.154802 -0.951018 -0.267584 -0.241356 -0.441811 -0.864031 -0.5 -0.309017 -0.809017 -0.43296 -0.587762 -0.683434 -0.714865 -0.149166 -0.683167 -0.683434 -0.43296 -0.587762 -0.587762 -0.683434 -0.43296 -0.864031 -0.241356 -0.441811 -0.809017 -0.5 -0.309017 -0.683167 -0.714865 -0.149166 -0.951018 -0.267584 -0.154802 -0.951018 -0.267584 0.154802 -0.855345 -0.518058 5.04833e-09 -0.864031 -0.241356 0.441811 -0.809017 -0.5 0.309017 -0.683167 -0.714865 0.149166 -0.714865 -0.149166 0.683167 -0.5 -0.309017 0.809017 -0.683434 -0.43296 0.587762 -0.241356 -0.441811 0.864031 -0.43296 -0.587762 0.683434 -0.587762 -0.683434 0.43296 -5.03782e-09 -0.292645 0.956221 -0.267584 -0.154802 0.951018 2.52357e-09 0 1 -0.518058 5.04833e-09 0.855345 -0.267584 0.154802 0.951018 0 0.292645 0.956221 0.518058 2.52416e-09 0.855345 0.267584 -0.154802 0.951018 0.267584 0.154802 0.951018 0 0.292645 -0.956221 2.52357e-09 0 -1 0.267584 0.154802 -0.951018 -2.51891e-09 -0.292645 -0.956221 0.267584 -0.154802 -0.951018 0.518058 -5.04833e-09 -0.855345 -0.267584 0.154802 -0.951018 -0.518058 -5.04833e-09 -0.855345 -0.267584 -0.154802 -0.951018</float_array>
                    <technique_common>
                        <accessor count="162" source="#ID7" stride="3">
                            <param name="X" type="float"/>
                            <param name="Y" type="float"/>
                            <param name="Z" type="float"/>
                        </accessor>
                    </technique_common>
                </source>
                <source id="ID8">
                    <float_array id="ID9" count="2">0 1</float_array>
                    <technique_common>
                        <accessor count="1" source="#ID9" stride="2">
                            <param name="S" type="float"/>
                            <param name="T" type="float"/>
                        </accessor>
                    </technique_common>
                </source>
                <vertices id="ID10">
                    <input semantic="POSITION" source="#ID4"/>
                </vertices>
                <triangles count="320" material="">
                    <input offset="0" semantic="VERTEX" source="#ID10"/>
                    <input offset="1" semantic="NORMAL" source="#ID6"/>
                    <input offset="2" semantic="TEXCOORD" source="#ID8" set="0"/>
                    <p>20 2 0 14 1 0 9 0 0 102 4 0 13 3 0 14 1 0 14 1 0 20 2 0 102 4 0 19 5 0 102 4 0 20 2 0 103 7 0 12 6 0 13 3 0 13 3 0 102 4 0 103 7 0 104 8 0 103 7 0 102 4 0 102 4 0 19 5 0 104 8 0 18 9 0 104 8 0 19 5 0 17 11 0 8 10 0 12 6 0 12 6 0 103 7 0 17 11 0 16 12 0 17 11 0 103 7 0 103 7 0 104 8 0 16 12 0 15 13 0 16 12 0 104 8 0 104 8 0 18 9 0 15 13 0 0 14 0 15 13 0 18 9 0 17 11 0 23 15 0 8 10 0 105 17 0 22 16 0 23 15 0 23 15 0 17 11 0 105 17 0 16 12 0 105 17 0 17 11 0 106 19 0 21 18 0 22 16 0 22 16 0 105 17 0 106 19 0 107 20 0 106 19 0 105 17 0 105 17 0 16 12 0 107 20 0 15 13 0 107 20 0 16 12 0 26 22 0 2 21 0 21 18 0 21 18 0 106 19 0 26 22 0 25 23 0 26 22 0 106 19 0 106 19 0 107 20 0 25 23 0 24 24 0 25 23 0 107 20 0 107 20 0 15 13 0 24 24 0 0 14 0 24 24 0 15 13 0 23 15 0 29 25 0 8 10 0 108 27 0 28 26 0 29 25 0 29 25 0 23 15 0 108 27 0 22 16 0 108 27 0 23 15 0 109 29 0 27 28 0 28 26 0 28 26 0 108 27 0 109 29 0 110 30 0 109 29 0 108 27 0 108 27 0 22 16 0 110 30 0 21 18 0 110 30 0 22 16 0 32 32 0 3 31 0 27 28 0 27 28 0 109 29 0 32 32 0 31 33 0 32 32 0 109 29 0 109 29 0 110 30 0 31 33 0 30 34 0 31 33 0 110 30 0 110 30 0 21 18 0 30 34 0 2 21 0 30 34 0 21 18 0 29 25 0 35 35 0 8 10 0 111 37 0 34 36 0 35 35 0 35 35 0 29 25 0 111 37 0 28 26 0 111 37 0 29 25 0 112 39 0 33 38 0 34 36 0 34 36 0 111 37 0 112 39 0 113 40 0 112 39 0 111 37 0 111 37 0 28 26 0 113 40 0 27 28 0 113 40 0 28 26 0 38 42 0 4 41 0 33 38 0 33 38 0 112 39 0 38 42 0 37 43 0 38 42 0 112 39 0 112 39 0 113 40 0 37 43 0 36 44 0 37 43 0 113 40 0 113 40 0 27 28 0 36 44 0 3 31 0 36 44 0 27 28 0 35 35 0 12 6 0 8 10 0 114 45 0 13 3 0 12 6 0 12 6 0 35 35 0 114 45 0 34 36 0 114 45 0 35 35 0 115 46 0 14 1 0 13 3 0 13 3 0 114 45 0 115 46 0 116 47 0 115 46 0 114 45 0 114 45 0 34 36 0 116 47 0 33 38 0 116 47 0 34 36 0 41 48 0 9 0 0 14 1 0 14 1 0 115 46 0 41 48 0 40 49 0 41 48 0 115 46 0 115 46 0 116 47 0 40 49 0 39 50 0 40 49 0 116 47 0 116 47 0 33 38 0 39 50 0 4 41 0 39 50 0 33 38 0 45 51 0 39 50 0 4 41 0 117 52 0 40 49 0 39 50 0 39 50 0 45 51 0 117 52 0 46 53 0 117 52 0 45 51 0 118 54 0 41 48 0 40 49 0 40 49 0 117 52 0 118 54 0 119 55 0 118 54 0 117 52 0 117 52 0 46 53 0 119 55 0 47 56 0 119 55 0 46 53 0 44 57 0 9 0 0 41 48 0 41 48 0 118 54 0 44 57 0 43 58 0 44 57 0 118 54 0 118 54 0 119 55 0 43 58 0 42 59 0 43 58 0 119 55 0 119 55 0 47 56 0 42 59 0 6 60 0 42 59 0 47 56 0 51 61 0 42 59 0 6 60 0 120 62 0 43 58 0 42 59 0 42 59 0 51 61 0 120 62 0 52 63 0 120 62 0 51 61 0 121 64 0 44 57 0 43 58 0 43 58 0 120 62 0 121 64 0 122 65 0 121 64 0 120 62 0 120 62 0 52 63 0 122 65 0 53 66 0 122 65 0 52 63 0 50 67 0 9 0 0 44 57 0 44 57 0 121 64 0 50 67 0 49 68 0 50 67 0 121 64 0 121 64 0 122 65 0 49 68 0 48 69 0 49 68 0 122 65 0 122 65 0 53 66 0 48 69 0 7 70 0 48 69 0 53 66 0 50 67 0 20 2 0 9 0 0 123 71 0 19 5 0 20 2 0 20 2 0 50 67 0 123 71 0 49 68 0 123 71 0 50 67 0 124 72 0 18 9 0 19 5 0 19 5 0 123 71 0 124 72 0 125 73 0 124 72 0 123 71 0 123 71 0 49 68 0 125 73 0 48 69 0 125 73 0 49 68 0 54 74 0 0 14 0 18 9 0 18 9 0 124 72 0 54 74 0 55 75 0 54 74 0 124 72 0 124 72 0 125 73 0 55 75 0 56 76 0 55 75 0 125 73 0 125 73 0 48 69 0 56 76 0 7 70 0 56 76 0 48 69 0 65 79 0 59 78 0 11 77 0 126 81 0 58 80 0 59 78 0 59 78 0 65 79 0 126 81 0 64 82 0 126 81 0 65 79 0 127 84 0 57 83 0 58 80 0 58 80 0 126 81 0 127 84 0 128 85 0 127 84 0 126 81 0 126 81 0 64 82 0 128 85 0 63 86 0 128 85 0 64 82 0 62 88 0 10 87 0 57 83 0 57 83 0 127 84 0 62 88 0 61 89 0 62 88 0 127 84 0 127 84 0 128 85 0 61 89 0 60 90 0 61 89 0 128 85 0 128 85 0 63 86 0 60 90 0 1 91 0 60 90 0 63 86 0 71 93 0 66 92 0 2 21 0 129 95 0 67 94 0 66 92 0 66 92 0 71 93 0 129 95 0 70 96 0 129 95 0 71 93 0 130 98 0 68 97 0 67 94 0 67 94 0 129 95 0 130 98 0 131 99 0 130 98 0 129 95 0 129 95 0 70 96 0 131 99 0 69 100 0 131 99 0 70 96 0 65 79 0 11 77 0 68 97 0 68 97 0 130 98 0 65 79 0 64 82 0 65 79 0 130 98 0 130 98 0 131 99 0 64 82 0 63 86 0 64 82 0 131 99 0 131 99 0 69 100 0 63 86 0 1 91 0 63 86 0 69 100 0 66 92 0 30 34 0 2 21 0 132 101 0 31 33 0 30 34 0 30 34 0 66 92 0 132 101 0 67 94 0 132 101 0 66 92 0 133 102 0 32 32 0 31 33 0 31 33 0 132 101 0 133 102 0 134 103 0 133 102 0 132 101 0 132 101 0 67 94 0 134 103 0 68 97 0 134 103 0 67 94 0 72 104 0 3 31 0 32 32 0 32 32 0 133 102 0 72 104 0 73 105 0 72 104 0 133 102 0 133 102 0 134 103 0 73 105 0 74 106 0 73 105 0 134 103 0 134 103 0 68 97 0 74 106 0 11 77 0 74 106 0 68 97 0 72 104 0 75 107 0 3 31 0 135 109 0 76 108 0 75 107 0 75 107 0 72 104 0 135 109 0 73 105 0 135 109 0 72 104 0 136 111 0 77 110 0 76 108 0 76 108 0 135 109 0 136 111 0 137 112 0 136 111 0 135 109 0 135 109 0 73 105 0 137 112 0 74 106 0 137 112 0 73 105 0 78 114 0 5 113 0 77 110 0 77 110 0 136 111 0 78 114 0 79 115 0 78 114 0 136 111 0 136 111 0 137 112 0 79 115 0 80 116 0 79 115 0 137 112 0 137 112 0 74 106 0 80 116 0 11 77 0 80 116 0 74 106 0 78 114 0 81 117 0 5 113 0 138 119 0 82 118 0 81 117 0 81 117 0 78 114 0 138 119 0 79 115 0 138 119 0 78 114 0 139 121 0 83 120 0 82 118 0 82 118 0 138 119 0 139 121 0 140 122 0 139 121 0 138 119 0 138 119 0 79 115 0 140 122 0 80 116 0 140 122 0 79 115 0 57 83 0 10 87 0 83 120 0 83 120 0 139 121 0 57 83 0 58 80 0 57 83 0 139 121 0 139 121 0 140 122 0 58 80 0 59 78 0 58 80 0 140 122 0 140 122 0 80 116 0 59 78 0 11 77 0 59 78 0 80 116 0 81 117 0 84 123 0 5 113 0 141 125 0 85 124 0 84 123 0 84 123 0 81 117 0 141 125 0 82 118 0 141 125 0 81 117 0 142 127 0 86 126 0 85 124 0 85 124 0 141 125 0 142 127 0 143 128 0 142 127 0 141 125 0 141 125 0 82 118 0 143 128 0 83 120 0 143 128 0 82 118 0 87 129 0 6 60 0 86 126 0 86 126 0 142 127 0 87 129 0 88 130 0 87 129 0 142 127 0 142 127 0 143 128 0 88 130 0 89 131 0 88 130 0 143 128 0 143 128 0 83 120 0 89 131 0 10 87 0 89 131 0 83 120 0 87 129 0 51 61 0 6 60 0 144 132 0 52 63 0 51 61 0 51 61 0 87 129 0 144 132 0 88 130 0 144 132 0 87 129 0 145 133 0 53 66 0 52 63 0 52 63 0 144 132 0 145 133 0 146 134 0 145 133 0 144 132 0 144 132 0 88 130 0 146 134 0 89 131 0 146 134 0 88 130 0 90 135 0 7 70 0 53 66 0 53 66 0 145 133 0 90 135 0 91 136 0 90 135 0 145 133 0 145 133 0 146 134 0 91 136 0 92 137 0 91 136 0 146 134 0 146 134 0 89 131 0 92 137 0 10 87 0 92 137 0 89 131 0 90 135 0 95 138 0 7 70 0 147 140 0 94 139 0 95 138 0 95 138 0 90 135 0 147 140 0 91 136 0 147 140 0 90 135 0 148 142 0 93 141 0 94 139 0 94 139 0 147 140 0 148 142 0 149 143 0 148 142 0 147 140 0 147 140 0 91 136 0 149 143 0 92 137 0 149 143 0 91 136 0 60 90 0 1 91 0 93 141 0 93 141 0 148 142 0 60 90 0 61 89 0 60 90 0 148 142 0 148 142 0 149 143 0 61 89 0 62 88 0 61 89 0 149 143 0 149 143 0 92 137 0 62 88 0 10 87 0 62 88 0 92 137 0 98 144 0 93 141 0 1 91 0 150 145 0 94 139 0 93 141 0 93 141 0 98 144 0 150 145 0 97 146 0 150 145 0 98 144 0 151 147 0 95 138 0 94 139 0 94 139 0 150 145 0 151 147 0 152 148 0 151 147 0 150 145 0 150 145 0 97 146 0 152 148 0 96 149 0 152 148 0 97 146 0 56 76 0 7 70 0 95 138 0 95 138 0 151 147 0 56 76 0 55 75 0 56 76 0 151 147 0 151 147 0 152 148 0 55 75 0 54 74 0 55 75 0 152 148 0 152 148 0 96 149 0 54 74 0 0 14 0 54 74 0 96 149 0 26 22 0 71 93 0 2 21 0 153 150 0 70 96 0 71 93 0 71 93 0 26 22 0 153 150 0 25 23 0 153 150 0 26 22 0 154 151 0 69 100 0 70 96 0 70 96 0 153 150 0 154 151 0 155 152 0 154 151 0 153 150 0 153 150 0 25 23 0 155 152 0 24 24 0 155 152 0 25 23 0 98 144 0 1 91 0 69 100 0 69 100 0 154 151 0 98 144 0 97 146 0 98 144 0 154 151 0 154 151 0 155 152 0 97 146 0 96 149 0 97 146 0 155 152 0 155 152 0 24 24 0 96 149 0 0 14 0 96 149 0 24 24 0 38 42 0 99 153 0 4 41 0 156 155 0 100 154 0 99 153 0 99 153 0 38 42 0 156 155 0 37 43 0 156 155 0 38 42 0 157 157 0 101 156 0 100 154 0 100 154 0 156 155 0 157 157 0 158 158 0 157 157 0 156 155 0 156 155 0 37 43 0 158 158 0 36 44 0 158 158 0 37 43 0 77 110 0 5 113 0 101 156 0 101 156 0 157 157 0 77 110 0 76 108 0 77 110 0 157 157 0 157 157 0 158 158 0 76 108 0 75 107 0 76 108 0 158 158 0 158 158 0 36 44 0 75 107 0 3 31 0 75 107 0 36 44 0 99 153 0 45 51 0 4 41 0 159 159 0 46 53 0 45 51 0 45 51 0 99 153 0 159 159 0 100 154 0 159 159 0 99 153 0 160 160 0 47 56 0 46 53 0 46 53 0 159 159 0 160 160 0 161 161 0 160 160 0 159 159 0 159 159 0 100 154 0 161 161 0 101 156 0 161 161 0 100 154 0 86 126 0 6 60 0 47 56 0 47 56 0 160 160 0 86 126 0 85 124 0 86 126 0 160 160 0 160 160 0 161 161 0 85 124 0 84 123 0 85 124 0 161 161 0 161 161 0 101 156 0 84 123 0 5 113 0 84 123 0 101 156 0</p>
                </triangles>
            </mesh>
        </geometry>
    </library_geometries>
    <library_visual_scenes>
        <visual_scene id="ID1">
            <node id="ID2" name="''' + name + '''">
                <translate sid="translate">0 0 -0</translate>
                <rotate sid="rotateY">0 1 0 -0</rotate>
                <rotate sid="rotateX">1 0 0 0</rotate>
                <rotate sid="rotateZ">0 0 1 -0</rotate>
                <scale sid="scale">1 1 1</scale>
                <instance_geometry url="#ID3"/>
            </node>
        </visual_scene>
    </library_visual_scenes>
    <scene>
        <instance_visual_scene url="#ID1"/>
    </scene>
</COLLADA>
''')

# pl = PDBList(pdb=pdbpath)
parser = PDBParser(PERMISSIVE=True, QUIET=True)  # QUIET=True suppresses warnings about PDB files


# if not os.path.isdir('dejavus'):
#    print 'making dejavus directory'
#    os.mkdir('dejavus')
#
# dejavupath = model_dir + 'dejavus' + os.sep


def writeIngredient(name, pdb, molarity, mw, color):
    print('write_Ingredient ' + name)
    pdbfn = pdbpath + str(pdb) + '.pdb'
    print(pdbfn)
    print(pdb)

    proxy = buildProxy(name, pdb, tree_sphere_radius, pdbfn, mw, surface=False, overwrite=True)

    positions = proxy[0]
    #    saveDejaVuMesh(name + '_cl', "", positions)      # don't think we need to save DejaVu file for clusters - that info is contained in the ingredient files.

    ingredient = open(model_dir + "%s.json" % name, "w")

    ingredient.write(str('{\n'))
    ingredient.write(str('    "packingMode": "random",\n'))
    ingredient.write(str('    "partners_position": [],\n'))
    ingredient.write(str('    "weight": 0.20000000000000001,\n'))
    ingredient.write(str('    "color": [' + color + '],\n'))
    # ingredient.write(str('        ' + str(random.random()) + ',\n'))  # random color
    # ingredient.write(str('        ' + str(random.random()) + ',\n'))
    # ingredient.write(str('        ' + str(random.random()) + '\n'))
    # ingredient.write(str('    ],\n'))
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
    ingredient.write(str('        0.1,\n'))
    ingredient.write(str('        0.1,\n'))  # what should jitterMax settings be?
    ingredient.write(str('        0.1\n'))
    ingredient.write(str('    ],\n'))
    ingredient.write(str('    "packingPriority": 0,\n'))
    ingredient.write(str('    "rotAxis": [\n'))
    ingredient.write(str('        0.0,\n'))
    ingredient.write(str('        2.0,\n'))
    ingredient.write(str('        1.0\n'))
    ingredient.write(str('    ],\n'))
    ingredient.write(str('    "overwrite_nbMol_value": 0,\n'))
    ingredient.write(str('    "nbJitter": 5,\n'))
    ingredient.write(str('    "molarity": ' + str(molarity) + ',\n'))
    ingredient.write(str('    "useOrientBias": false,\n'))
    ingredient.write(str('    "rotRange": 6.2831000000000001,\n'))
    ingredient.write(str('    "coordsystem": "left",\n'))
    ingredient.write(str('    "meshFile": "' + name + "_coarse.dae" + '",\n'))
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
    ingredient.write(str('        ' + str(positions) + '\n'))
    ingredient.write(str('    ],\n'))
    ingredient.write(str('    "meshName": "' + name + '",\n'))
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
    recipe.write(str('  "smallestProteinSize":10,\n'))
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
    recipe.write(str(' },'))

    recipe.write(str('''
 "cytoplasme":{
  "ingredients":{
  }
 },
 "compartments":{\n'''))

    compartments = []
    for header in headers:
        if header[0] == 'C' and header[1] == '_':
            compartments.append(header)
            print(header + ' compartment detected')
    print('compartments total:')
    print(compartments)

    first_compartment = True

    for compartment in compartments:

        if not first_compartment:  # the first won't be preceeded by a comma. (the last can't have a comma)
            recipe.write(str(',\n'))
        first_compartment = False

        recipe.write(str('''  "''' + str(compartment[2:]) + '''":{
   "geom":"''' + str(compartment[2:]) + '''.dae",
   "name":"''' + str(compartment[2:]) + '''",
   "surface":{
   "ingredients":{}
   },
   "interior":{
    "ingredients":{
'''))

        first = True

        for x in range(1, len(all_data)):  # for each entry in the input .csv
            if all_data[x][headers['INCLUDE']]:
                mw = ''  # assuming average human protein size of 320 aa, 110 Da per aa
                name = ''
                pdb = ''
                color = ''
                if 'NAME' in headers:
                    name = str(all_data[x][headers['NAME']])
                if not name:
                    name = 'Unnamed_ingredient_1'
                if 'PDB' in headers:
                    pdb = all_data[x][headers['PDB']]
                if not pdb:
                    pdb = 'null'
                if 'MW' in headers:
                    mw = all_data[x][headers['MW']]
                if not mw:
                    mw = '35200'  # assuming size of 320 aa
                mw = int(mw)
                # if 'MOLARITY' in headers:
                #     molarity = all_data[x][headers['MOLARITY']]
                # if not molarity:
                #     molarity = 0
                if 'COLOR' in headers:
                    color = all_data[x][headers['COLOR']]
                if not color:
                    color = str(random.random()) + ', ' + str(random.random()) + ', ' + str(random.random())
                molarity = all_data[x][headers[compartment]]
                if not molarity:
                    continue

                # pdb = pdb.split("_")[0]  # gets rid of any "_X" after PDB name
                # if not pdb:
                #     pdb = 'pdb' + str(x)
                print('pdb = ' + pdb)
                if not pdb == 'null' and not os.path.isfile(pdbpath + str(pdb) + '.pdb'):
                    write_pdbFile(pdb, pdbpath)  # download it to PDB directory
                if not pdb == 'null' and not os.path.isdir(model_dir + 'PDB/' + pdb):  # if the PDB's specific directory is not already there
                    os.mkdir(model_dir + 'PDB/' + pdb)
                    print('making directory for ' + pdb)
                    if not os.path.isfile(pdbpath + str(pdb) + os.sep + str(pdb) + '.pdb'):
                        write_pdbFile(pdb, str(pdbpath + pdb + os.sep))  # download it to its directory
                    # if not newfile:
                    #     print(pdb + ' NOT FOUND IN PDB')
                    #     print('NO PDB OR MW - making 320 aa dummy sphere for ' + pdb)
                    #     continue
                name = nameFix(name)
                if not first:  # the first won't be preceeded by a comma. (the last can't have a comma)
                    recipe.write(str(',\n'))
                first = False
                recipe.write(str('     "' + name + '":{\n'))
                recipe.write(str('      "include":"' + name + '.json",\n'))
                recipe.write(str('      "name":"' + name + '"\n'))
                recipe.write(str('     }'))

                writeIngredient(name, pdb, molarity, mw, color)

        recipe.write(str('\n    }\n'))  # end of ingredients for interior of this compartment
        recipe.write(str('   }\n'))  # end of interior of this compartment

        recipe.write(str('  }'))  # end of this compartment

    recipe.write(str('\n }\n'))  # end of "compartments"
    recipe.write(str('}'))  # end of object

    recipe.close()


writeRecipe()

os.system('say "Beer time."')

print("done")
