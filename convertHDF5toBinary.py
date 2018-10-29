import numpy as np
import struct
import sys
import os
import h5py

"""
------------------------------------------

Execute by typing:
python convertHDF5toBinary.py snaplist.txt

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
input files, but without the .hdf5 extension. You can easily change that
in the MAIN part on the bottom of this file.

The format-2 option might not be fully working, it seems to be working, 
but I'm not sure how to test it.

The catalogs and header are the default from GADGET2. If you want anything
extra you can add it in the arrays below.

"""

#-----------Adjust this part if necessary-------------------------------------------------------------

#If there are any extra HDF5 datasets you'd like to convert (or if they have different names), you can add them here:
keysHDF5 = (['Coordinates', 'Velocities', 'ParticleIDs', 'Masses', 'InternalEnergy', 'Density', 'SmoothingLength', 
  'Potential', 'Acceleration', 'RateOfChangeOfEntropy', 'TimeStep'])
#Their binary names are suitable for converting to format-2, the order should be the same as keysHDF5 and 
#is written out in binary in this order:
keysBin = (['POS ', 'VEL ', 'ID  ', 'MASS', 'U   ', 'RHO ', 'HSML', 'POT ', 'ACCE', 'ENDT', 'TSTP'])
#The binary c-types the data is being written out to:
datatype = ['f', 'f', 'i', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f']
#The HDF5 header names, you can add stuff if you need to:
keysHDF5Header = (['NumPart_ThisFile', 'MassTable', 'Time', 'Redshift', 'Flag_Sfr', 'Flag_Feedback', 'NumPart_Total', 
  'Flag_Cooling', 'NumFilesPerSnapshot', 'BoxSize', 'Omega0', 'OmegaLambda', 'HubbleParam', 'Flag_StellarAge', 'Flag_Metals', 
  'NumPart_Total_HighWord'])
#The binary c-types of the header attributes
datatypeBinHeader = "".join(['i'*6, 'd'*6, 'd', 'd', 'i', 'i', 'i'*6, 'i', 'i', 'd', 'd', 'd', 'd', 'i', 'i', 'i'*6])

#------------------------------------------------------------------------------------------------------


keysmatch = {}
datatypes = {}
for i in range(len(keysHDF5)):
  keysmatch[keysHDF5[i]] = keysBin[i]
  datatypes[keysBin[i]] = datatype[i]



def read_hdf5file(snapfile):
  d_dataset = {}
  d_header = {}

  #Opening file
  if os.path.isfile(snapfile):
    dataopen = h5py.File(snapfile, "r")
  else:
    sys.exit("Error: Could not open %s" %(snapfile))
  
  #Reading header
  Header = dataopen['Header']
  for key in Header.attrs.keys():
    d_header[key] = Header.attrs[key]
  # d_header['MassTable'][0] = 0.0
  # d_header['NumPart_ThisFile'][0] = 0
  # d_header['NumPart_Total'][0] = 0
  # d_header['NumPart_Total_HighWord'][0] = 0
  
  #Reading datasets
  for i in range(6):
    #Checking if the particle type exists
    particle_type = 'PartType{}'.format(int(i)) in dataopen
    data_sets = []

    #If it does, a data set will be made and written
    if particle_type:
      d_dataset[i] = {}
      for j in dataopen['PartType{}'.format(int(i))].id:
        d_dataset[i][keysmatch[j.decode('utf-8')]] = dataopen['PartType{}/'.format(int(i)) + j.decode('utf-8')][:]

  return d_header, d_dataset
  

def pack_header(d_header):
  h_data = []
  for key in keysHDF5Header:
    if isinstance(d_header[key], (list, tuple, np.ndarray)):
      for i in range(len(d_header[key])):
        h_data.append(d_header[key][i])
    else:
      h_data.append(d_header[key])

  # Fill up to 256 bytes
  for i in np.arange(256 - struct.calcsize(datatypeBinHeader)):
    h_data.append('\0'.encode('utf-8'))
  s = struct.Struct(datatypeBinHeader + 'c'*(256-struct.calcsize(datatypeBinHeader)))
  packed_data = s.pack(*h_data)
  return packed_data


def write_dummy(f, values_list):
  for i in values_list:
    dummy = [i]
    s = struct.Struct('i')
    d = s.pack(*dummy)
    f.write(d)


def write_file(filename, d_header, d_dataset, format2=False):
  uit = open(filename, "wb")

  if format2:
    write_dummy(uit, [8])
    uit.write(struct.pack('c' * 4, *[i.encode('utf-8') for i in 'HEAD']))
    write_dummy(uit, [8])

  nbytes = 256
  write_dummy(uit, [nbytes]) 
  uit.write(pack_header(d_header))
  write_dummy(uit, [nbytes])

  #Loop over all given keys
  for key in keysBin:
    nbytes = 0
    #See how many bytes have to be allocated in dummy
    for parttype in d_dataset.keys():
      if key in d_dataset[parttype].keys():
        if isinstance(d_dataset[parttype][key][0], (list, tuple, np.ndarray)):
          nbytes += len(d_dataset[parttype][key]) * struct.calcsize(datatypes[key]) * 3
        else:
          nbytes += len(d_dataset[parttype][key]) * struct.calcsize(datatypes[key])
      else:
        continue

    #Check if the block exists, otherwise skip
    if nbytes == 0:
      continue

    if format2:
      write_dummy(uit, [8])
      uit.write(struct.pack('c'*4, *[i.encode('utf-8') for i in key]))
      write_dummy(uit, [8])

    #Write size of the dummy
    print(nbytes)
    write_dummy(uit, [nbytes])

    #Write particle data
    for parttype in d_dataset.keys():
      if isinstance(d_dataset[parttype][key][0], (list, tuple, np.ndarray)):
        uit.write(struct.pack(datatypes[key]*3*len(d_dataset[parttype][key]), *d_dataset[parttype][key][:].flatten()))
      else:
        uit.write(struct.pack(datatypes[key]*len(d_dataset[parttype][key]), *d_dataset[parttype][key][:]))
    #Write buffer dummy
    write_dummy(uit, [nbytes])
    print('Finished writing dataset: ', key)

  uit.close()


def convert_file(infile, outfile, format2=False):
  d_header, d_dataset = read_hdf5file(infile)
  write_file(outfile, d_header, d_dataset)


#--------------MAIN--------------------------------------------------------------
filenames = np.genfromtxt(sys.argv[1], dtype='str')

indir = filenames[0]
outdir = filenames[1]

for i in range(2, len(filenames)):
  print('convert ', indir+filenames[i], ' to ', outdir+filenames[i][:-5])
  convert_file(indir+filenames[i], outdir+filenames[i][:-5])

