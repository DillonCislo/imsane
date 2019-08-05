import h5py
import argparse
import numpy as np

"""Unpack a N-Dim hdf5 file into an (N-1)-Dim HDF5 file.

Example Usage
-------------
python unstack_hdf5.py -fn my_dataset.h5 -dset datasetname -outdset new_dataset -out my_dataset_wildcard.h5 -indexsz 5

"""

parser = argparse.ArgumentParser(description='Unpack an HDF5 file into separate HDF5 files on disk.')

# Options
parser.add_argument('-fn', '--data_filename', help='The filename to load', type=str, default='empty_string')
parser.add_argument('-dset', '--dataset_name', help='The name of the dataset inside the HDF5 file',
                    type=str, default='exported_data')
parser.add_argument('-outdset', '--output_dataset_name',
                    help='The name of the output dataset inside the output HDF5 file',
                    type=str, default='exported_data')
parser.add_argument('-dim', '--unpacking_dimension', help='The dimension of the dataset to unpack', type=int, default=0)
parser.add_argument('-out', '--output_filename',
                    help="The filename base to write data to on disk, with 'wildcard' in place of the index",
                    type=str, default='empty_string')
parser.add_argument('-indexsz', '--indexsz', help='How many digits to use for the timepoint index (with leading zeros)',
                    type=int, default=5)

args = parser.parse_args()

if args.data_filename in ['empty_string', 'none', '']:
    raise RuntimeError('Please supply an HDF5 filename to unpack')


# First load the dataset from the HDF5 file
f = h5py.File(args.data_filename, 'r')
dset = f[args.dataset_name][:]
f.close()
print('np.shape(dset) = ', np.shape(dset))

# Unpack the dataset
slc = [slice(None)] * len(np.shape(dset))
slc[args.unpacking_dimension] = 0
lengthdset = len(np.shape(dset))

for ii in range(np.shape(dset)[args.unpacking_dimension]):
    print('ii = ' + str(ii))
    # First make a list for slicing all the axes
    slc = [slice(None)] * lengthdset
    slc[args.unpacking_dimension] = ii

    # Take this 'frame' from the stack
    # Could use numpy's take, but this copies things
    # frame = dset.take(indices=ii, axis=args.unpacking_dimension)
    frame = dset[tuple(slc)]

    # Save this (N-1)D frame to disk
    if 'wildcard' in args.output_filename:
        indexstr = '{0:0d}'.format(ii).rjust(args.indexsz, '0')
        outfn = args.output_filename.replace('wildcard', indexstr)
    else:
        outfn = args.output_filename
        print('WARNING: Potential overwrite problem: writing to ' + outfn)

    f = h5py.File(outfn, 'w')
    outdset = f.create_dataset(args.output_dataset_name, data=frame, chunks=True)
    f.close()



