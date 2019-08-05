import h5py
import argparse


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

# Unpack the dataset
slc = [slice(None)] * len(dset.shape)
slc[args.unpacking_dimension] = 0

for ii in range(len(dset[slc])):
    # First make a list for slicing all the axes
    slc = [slice(None)] * len(dset.shape)
    slc[args.unpacking_dimension] = ii
    # Take this 'frame' from the stack
    frame = dset[slc]

    # Save this (N-1)D frame to disk
    if 'wildcard' in args.output_filename:
        indexstr = print(ii.zfill(args.indexsz))
        outfn = args.output_filename.replace('wildcard', indexstr)
    else:
        outfn = args.output_filename
        print('WARNING: Potential overwrite problem: writing to ' + outfn)

    f = h5py.File(args.outfn, 'w')
    dset = f.create_dataset(args.output_dataset_name, data=frame, chunks=True)
    f.close()



