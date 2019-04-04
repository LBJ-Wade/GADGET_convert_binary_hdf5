import numpy as np
import sys
import os
import h5py

"""
------------------------------------------

Execute by typing:
python convertBinarytoHDF5.py snaplist.txt

------------------------------------------
The snaplist.txt should look like this:
/path/to/input/directory/
/path/to/output/directory/
snapname_1.hdf5
snapname_2.hdf5
.
.
.
snapname_N.hdf5

-----------Notes--------------------------------------------------------
I've set it up so that the outputfiles have the same name as the 
input files, but now with a .hdf5 extension. You can easily change that
in the MAIN part on the bottom of this file.

The format-2 option is not working, but it should if you increase the 'dummylength'
in read_binaryfile.

The catalogs and header are the default from GADGET2. If you want anything
extra you can add it to the arrays below.

In principle there is no need to remove keys from the lists below. Properties are read in
using the header information from the binary file.
All particle types will read in Coordinates, Velocities, IDs, Masses in case of 
its entry in MassTable being zero.

PartType0 will read the InternalEnergy, Density and SmoothingLength in addition. If you
are converting ICs that lack any of this information, you can set the flags missingDensity,
missingSmoothingLength, missingInternalEnergy to True.

If dealing with long ints, you can set the longInts flag as True (I haven't tested this option yet).

"""

#-----------Adjust this part if necessary-------------------------------------------------------------
#Some flags (explained in 'Notes')
d_flags = {
"longInts" : False,
"missingDensity" : False,
"missingSmoothingLength" : False,
"missingInternalEnergy" : False
}

#If there are any extra HDF5 datasets you'd like to convert (or if you want to give them different names), you can add them here.
#If you need to change the order, you will need to change that in read_binaryfile too.
keysHDF5 = (['Coordinates', 'Velocities', 'ParticleIDs', 'Masses', 'InternalEnergy', 'Density', 'SmoothingLength', 
  'Potential', 'Acceleration', 'RateOfChangeOfEntropy', 'TimeStep'])

#The HDF5 header names, you can add stuff if you need to, but don't forget to add it too to the datatypeHeader and countsHeader
#below:
keysHDF5Header = (['NumPart_ThisFile', 'MassTable', 'Time', 'Redshift', 'Flag_Sfr', 'Flag_Feedback', 'NumPart_Total', 
  'Flag_Cooling', 'NumFilesPerSnapshot', 'BoxSize', 'Omega0', 'OmegaLambda', 'HubbleParam', 'Flag_StellarAge', 'Flag_Metals', 
  'NumPart_Total_HighWord', 'Flag_Entropy_ICs'])

if d_flags["longInts"]:
  datatypeHeader = [np.uint64, np.float64, np.float64, np.float64, np.int32, np.int32, np.uint64, np.int32, np.int32, np.float64, 
    np.float64, np.float64, np.float64, np.int32, np.int32, np.uint64, np.int32]
else:
  datatypeHeader = [np.uint32, np.float64, np.float64, np.float64, np.int32, np.int32, np.uint32, np.int32, np.int32, np.float64, 
    np.float64, np.float64, np.float64, np.int32, np.int32, np.uint32, np.int32]
countsHeader = [6, 6, 1, 1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1]

#------------------------------------------------------------------------------------------------------

#----------Not implemented yet-------------
# #Their binary names are suitable for converting to format-2, the order should be the same as keysHDF5 and 
# #is written out in binary in this order:
# keysBin = (['POS ', 'VEL ', 'ID  ', 'MASS', 'U   ', 'RHO ', 'HSML', 'POT ', 'ACCE', 'ENDT', 'TSTP'])
# keysmatch = {}
# for i in range(len(keysHDF5)):
#   keysmatch[keysBin[i]] = keysHDF5[i]


def createheader(dataopen, d_header):
  header = dataopen.create_group("Header")
  i = 0
  for key in d_header.keys():
    if countsHeader[i] > 1:
      header.attrs.__setitem__(key, d_header[key].astype(dtype=datatypeHeader[i]))
    else:
      header.attrs.__setitem__(key, datatypeHeader[i](d_header[key]))
    i += 1

def write_hdf5file(snapfile, d_header, d_dataset):
  print("Writing data to "+snapfile)
  dataopen = h5py.File(snapfile, 'w')
  createheader(dataopen, d_header)
  for parttype in d_dataset.keys():
    partset = dataopen.create_group("PartType{}".format(int(parttype)))
    for key in d_dataset[parttype].keys():
      partset.create_dataset(key, data = d_dataset[parttype][key])
  print("Finished!")

def read_binaryfile(filename, d_flags):
  d_dataset = {}

  print("Reading "+filename+"...")
  f = open(filename, "rb")

  dummylength = 4

  d_header = {}

  HeaderLength = np.fromfile(f,dtype=np.uint32,count=1)

  i = 0

  bn = 0
  for key in keysHDF5Header:
    bn += np.dtype(datatypeHeader[i]).itemsize*countsHeader[i]
    d_header[key] = np.fromfile(f, dtype = datatypeHeader[i], count = countsHeader[i])
    i += 1
  
  #Read the rest of the empty header
  f.read(int(HeaderLength - bn))

  #Print the header to screen
  print("----------------HEADER INFO---------------------")
  for key in d_header.keys():
    if len(d_header[key]) == 1:
      print(key, '\t', d_header[key][0])
    else:
      print(key, '\t', d_header[key])
  print("------------------------------------------------")

  f.read(dummylength)

  #Construct dataset skeleton based on header information
  print("Constructing data skeleton...")
  for parttype in range(6):
    if d_header['NumPart_ThisFile'][parttype] > 0:
      print('Particle type %i present' %parttype)
      d_dataset[parttype] = {}
      for i in range(2):
        d_dataset[parttype][keysHDF5[i]] = np.zeros((d_header['NumPart_ThisFile'][parttype], 3), dtype=np.float32)
      if d_flags['longInts']:
        d_dataset[parttype][keysHDF5[2]] = np.zeros(d_header['NumPart_ThisFile'][parttype], dtype=np.uint64)
      else:
        d_dataset[parttype][keysHDF5[2]] = np.zeros(d_header['NumPart_ThisFile'][parttype], dtype=np.uint32)

      #If no mass specified in the MassTable, it is read in per particle
      if d_header['MassTable'][parttype] == 0:
        d_dataset[parttype][keysHDF5[3]] = np.zeros(d_header['NumPart_ThisFile'][parttype], dtype=np.float32)
      
      #Gas particles have as a default internal energy, density and smoothing lengths.
      if parttype == 0:
        if d_flags['missingInternalEnergy'] == False:
          d_dataset[parttype][keysHDF5[4]] = np.zeros(d_header['NumPart_ThisFile'][parttype], dtype=np.float32)
        if d_flags['missingDensity'] == False:
          d_dataset[parttype][keysHDF5[5]] = np.zeros(d_header['NumPart_ThisFile'][parttype], dtype=np.float32)
        if d_flags['missingSmoothingLength'] == False:
          d_dataset[parttype][keysHDF5[6]] = np.zeros(d_header['NumPart_ThisFile'][parttype], dtype=np.float32)
        if d_header['Flag_Entropy_ICs']:
          d_dataset[parttype][keysHDF5[9]] = np.zeros(d_header['NumPart_ThisFile'][parttype], dtype=np.float32)
  print("Finished skeleton")
  
  #Read Coordinates and Velocities
  print("Reading coordinates and velocities")
  for key in keysHDF5[0:2]:
    f.read(dummylength)
    for parttype in d_dataset.keys():
      for i in range(d_header['NumPart_ThisFile'][parttype]):
        for j in range(3):
          d_dataset[parttype][key][i][j] = np.fromfile(f,dtype=np.float32,count=1)
    f.read(dummylength)

  #Read IDs
  print("Reading IDs")
  key=keysHDF5[2]
  f.read(dummylength)
  for parttype in d_dataset.keys():
    if key in d_dataset[parttype].keys():
      for i in range(d_header['NumPart_ThisFile'][parttype]):
        if d_flags['longInts']:
          d_dataset[parttype][key][i] = np.fromfile(f,dtype=np.uint64,count=1)
        else:
          d_dataset[parttype][key][i] = np.fromfile(f,dtype=np.uint32,count=1)
  f.read(dummylength)

  #Read Mass and gas properties
  for key in keysHDF5[3:]:
    dummydone = False
    for parttype in d_dataset.keys():
      if key in d_dataset[parttype].keys():
        print("Reading",key,"for particle type %i"%parttype)
        if not dummydone:
          f.read(dummylength)
          dummydone = True
        for i in range(d_header['NumPart_ThisFile'][parttype]):
          d_dataset[parttype][key][i] = np.fromfile(f,dtype=np.float32,count=1)
    if dummydone:
      f.read(dummylength)
  
  f.close()
  print("Finished reading binary file (/-_-)/")
  return d_header,d_dataset


def convert_file(infile, outfile, d_flags):
  d_header, d_dataset = read_binaryfile(infile, d_flags)
  write_hdf5file(outfile, d_header, d_dataset)


#--------------MAIN--------------------------------------------------------------
filenames = np.genfromtxt(sys.argv[1], dtype='str')

indir = filenames[0]
outdir = filenames[1]

for i in range(2, len(filenames)):
  print('convert ', indir+filenames[i], ' to ', outdir+filenames[i]+'.hdf5')
  convert_file(indir+filenames[i], outdir+filenames[i]+'.hdf5', d_flags)

