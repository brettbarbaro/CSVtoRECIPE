"""
Created on Saturday, July 16, 2016
@author: Jared Truong
"""
import urllib
import os

# given pdb ID
# returns pdb file information
def fetch_pdb(pdbid):
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid.upper()
    return urllib.urlopen(url).read()

# given pdb ID and path of folder to store files('C:\\Users\\User\\Desktop\\pdbFiles')
# writes pdb file into given location with file name "pdbid.pdb"
def write_pdbFile(pdbid, locationPath):
    file = open(str(locationPath) + '\\' + str(pdbid) + '.pdb', 'w')
    data = fetch_pdb(pdbid)
    file.write(data)
    file.close()

# given list of pdbid strings, and path of folder to store files
# writes a pdb file for each given IDs
def savePdbFiles(pdbidList, locationPath):
    for i in range(len(pdbidList)):
        write_pdbFile(pdbidList[i], locationPath)

# given path of folder that stores pdb files
# returns list of strings with file paths
def fetch_pdbFileList(pdbFileDirectoryPath):
    pdbFileList = os.listdir(pdbFileDirectoryPath)
    for i in range(len(pdbFileList)):
        pdbFileList[i] = pdbFileDirectoryPath + '\\' + pdbFileList[i]
    return pdbFileList

