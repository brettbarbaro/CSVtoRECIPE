# all information taken from Stanford websites on Mycoplasm genitalium
# http://www.wholecellsimdb.org/simulation_batch/1
# http://www.wholecellsimdb.org/simulation_batch/2
# http://www.wholecellkb.org
# goal to extract copy numbers for protein complexes from simulation database in the form of h5 files
# each simulation contains 201 protein complexes in six celluar compartments, six simulation states
# and approximately 30000 frames (to represent different parts of a given simulation)
import h5py
import numpy as np
import xlwt
from tempfile import TemporaryFile

# given h5 file location, frame, and ignore_null boolean (ignore_null declares whether or no to ignore data that is 0)
# file location as a string (i.e. 'C:\\Users\\User\\Desktop\\Project\\simulation.h5')
# maximum frame value varies by simulation and is about 30000
# returns 3D array with information on compartment, name, and value
def getRawData(filename, frame=0, ignore_null=False):
    data = {}
    f = h5py.File(filename)
    dataComplex = f.get("states/ProteinComplex/counts/data")
    labelComplex = f.get("states/ProteinComplex/counts/labels/0")
    for i in range(6):
        if ignore_null:
            names = labelComplex[dataComplex[:, i, frame] > 0]
            counts = dataComplex[:, i, frame][dataComplex[:, i, frame] > 0]
        else:
            names = labelComplex
            counts = dataComplex[:, i, frame]
        compl = np.column_stack([names, counts])
        data[i] = np.vstack([compl]).tolist()
    return data

# returns 2D array with data on cytosol compartment given h5 file location
def getCData(filename, frame=0, ignore_null=False):
    return getRawData(filename, frame, ignore_null)[0]

# returns 2D array with data on dna compartment given h5 file location
def getDData(filename, frame=0, ignore_null=False):
    return getRawData(filename, frame, ignore_null)[1]

# returns 2D array with data on extracellular space compartment given h5 file location
def getESData(filename, frame=0, ignore_null=False):
    return getRawData(filename, frame, ignore_null)[2]

# returns 2D array with data on membrane compartment given h5 file location
def getMData(filename, frame=0, ignore_null=False):
    return getRawData(filename, frame, ignore_null)[3]

# returns 2D array with data on terminal organelle ctysol compartment given h5 file location
def getTOCData(filename, frame=0, ignore_null=False):
    return getRawData(filename, frame, ignore_null)[4]

# returns 2D array with data on terminal organelle membrane compartment given h5 file location
def getTOMData(filename, frame=0, ignore_null=False):
    return getRawData(filename, frame, ignore_null)[5]

# there are six simulation states in this data
# NASCENT = 0
# MATURE = 1
# INACTIVATED = 2
# BOUND = 3
# MISFOLDED = 4
# DAMAGED = 5
# given h5file location, frame, ignore_null boolean value, simulation state (as listed above)
# returns 3D array with data on specific simulation state
def getSimulation(filename, frame=0, ignore_null=False, simulation=1):
    if simulation == 0:
        startIndex = 0
    elif simulation == 1:
        startIndex = 201
    elif simulation == 2:
        startIndex = 402
    elif simulation == 3:
        startIndex = 603
    elif simulation == 4:
        startIndex = 804
    elif simulation == 5:
        startIndex = 1005
    else:
        startIndex = -1
    if startIndex >= 0:
        retMatrix = [getRawData(filename, frame, ignore_null)[i][startIndex:startIndex+201][0:201] for i in range(0, 6)]
    else:
        retMatrix = ['simulation parameter not valid']
    return retMatrix

# there is 6 compartments in this data
# c Cytosol 0
# d DNA     1
# e Extracellular Space  2
# m Membrane  3
# tc Terminal Organelle Cytosol  4
# tm Terminal Organelle Membrane  5
# localisation = {"c": 0, "d": 1, "e": 2, "m": 3, "tc": 4, "tm": 5}
# given raw data in a 3D array (see getRawData(filename, frame=0, ignore_null=False) function)
# creates new xls workbook with six sheets representing different compartments
# each sheet contains complex name and it respective count value
def dataToXls(data, fileName):
    book = xlwt.Workbook()
    sheet1 = book.add_sheet('cytosolData', cell_overwrite_ok=True)
    sheet2 = book.add_sheet('dnaData', cell_overwrite_ok=True)
    sheet3 = book.add_sheet('extracellularSpaceData', cell_overwrite_ok=True)
    sheet4 = book.add_sheet('membraneData', cell_overwrite_ok=True)
    sheet5 = book.add_sheet('terminalOrganelleCytosolData', cell_overwrite_ok=True)
    sheet6 = book.add_sheet('terminalOrganelleMembraneData', cell_overwrite_ok=True)

    cytosolData = data[0]
    dnaData = data[1]
    extracellularSpaceData = data[2]
    membraneData = data[3]
    terminalOrganelleCytosolData = data[4]
    terminalOrganelleMembraneData = data[5]

    for row, list in enumerate(cytosolData):
        for col, value in enumerate(list):
            sheet1.write(row, col, value)

    for row, list in enumerate(dnaData):
        for col, value in enumerate(list):
            sheet2.write(row, col, value)

    for row, list in enumerate(extracellularSpaceData):
        for col, value in enumerate(list):
            sheet3.write(row, col, value)

    for row, list in enumerate(membraneData):
        for col, value in enumerate(list):
            sheet4.write(row, col, value)

    for row, list in enumerate(terminalOrganelleCytosolData):
        for col, value in enumerate(list):
            sheet5.write(row, col, value)

    for row, list in enumerate(terminalOrganelleMembraneData):
        for col, value in enumerate(list):
            sheet6.write(row, col, value)

    name = fileName+".xls"
    book.save(name)
    book.save(TemporaryFile())

# given array of file names to be added to xls file, total number of files, name of xls output, frame, ignore_null value
# generates new xls workbook with six sheets representing different compartments
# each sheet containing complex name in the first column and data from all inputted simulations in subsequent columns
def multipleJ5ToXls(listOfFiles, numberOfIterations, fileName, frame=0, ignore_null=False):
    cytosolData = []
    dnaData = []
    extracellularSpaceData = []
    membraneData = []
    terminalOrganelleCytosolData = []
    terminalOrganelleMembraneData = []
    # separates data into 3D arrays representing specific compartments
    for i in range(numberOfIterations):
        cytosolData.append(getRawData(listOfFiles[i], frame, ignore_null)[0])
        dnaData.append(getRawData(listOfFiles[i], frame, ignore_null)[1])
        extracellularSpaceData.append(getRawData(listOfFiles[i], frame, ignore_null)[2])
        membraneData.append(getRawData(listOfFiles[i], frame, ignore_null)[3])
        terminalOrganelleCytosolData.append(getRawData(listOfFiles[i], frame, ignore_null)[4])
        terminalOrganelleMembraneData.append(getRawData(listOfFiles[i], frame, ignore_null)[5])
    # adds six sheets to new book
    book = xlwt.Workbook()
    sheet1 = book.add_sheet('cytosolData', cell_overwrite_ok=True)
    sheet2 = book.add_sheet('dnaData', cell_overwrite_ok=True)
    sheet3 = book.add_sheet('extracellularSpaceData', cell_overwrite_ok=True)
    sheet4 = book.add_sheet('membraneData', cell_overwrite_ok=True)
    sheet5 = book.add_sheet('terminalOrganelleCytosolData', cell_overwrite_ok=True)
    sheet6 = book.add_sheet('terminalOrganelleMembraneData', cell_overwrite_ok=True)
    # for each sheet, adds values starting in the second column,
    # then adds protein complex name labels in the first column
    for x, plane in enumerate(cytosolData):
        for row, list in enumerate(plane):
            for col, value in enumerate(list):
                if col == 1:
                    sheet1.write(row, col * x + 1, value)
    for x, plane in enumerate(cytosolData):
        for row, list in enumerate(plane):
            for col, value in enumerate(list):
                if col == 0:
                    sheet1.write(row, col, value)
        break

    for x, plane in enumerate(dnaData):
        for row, list in enumerate(plane):
            for col, value in enumerate(list):
                if col == 1:
                    sheet2.write(row, col * x + 1, value)
    for x, plane in enumerate(dnaData):
        for row, list in enumerate(plane):
            for col, value in enumerate(list):
                if col == 0:
                    sheet2.write(row, col, value)
        break

    for x, plane in enumerate(extracellularSpaceData):
        for row, list in enumerate(plane):
            for col, value in enumerate(list):
                if col == 1:
                    sheet3.write(row, col * x + 1, value)
    for x, plane in enumerate(extracellularSpaceData):
        for row, list in enumerate(plane):
            for col, value in enumerate(list):
                if col == 0:
                    sheet3.write(row, col, value)
        break

    for x, plane in enumerate(membraneData):
        for row, list in enumerate(plane):
            for col, value in enumerate(list):
                if col == 1:
                    sheet4.write(row, col * x + 1, value)
    for x, plane in enumerate(membraneData):
        for row, list in enumerate(plane):
            for col, value in enumerate(list):
                if col == 0:
                    sheet4.write(row, col, value)
        break

    for x, plane in enumerate(terminalOrganelleCytosolData):
        for row, list in enumerate(plane):
            for col, value in enumerate(list):
                if col == 1:
                    sheet5.write(row, col * x + 1, value)
    for x, plane in enumerate(cytosolData):
        for row, list in enumerate(plane):
            for col, value in enumerate(list):
                if col == 0:
                    sheet5.write(row, col, value)
        break

    for x, plane in enumerate(terminalOrganelleMembraneData):
        for row, list in enumerate(plane):
            for col, value in enumerate(list):
                if col == 1:
                    sheet6.write(row, col * x + 1, value)
    for x, plane in enumerate(terminalOrganelleMembraneData):
        for row, list in enumerate(plane):
            for col, value in enumerate(list):
                if col == 0:
                    sheet6.write(row, col, value)
        break
    # saves file
    name = fileName+".xls"
    book.save(name)
    book.save(TemporaryFile())

# 20160713 work in progress...