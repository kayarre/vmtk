## Module:    $RCSfile: vmtknumpywriter.py,v $
## Language:  Python
## Date:      June 10, 2017
## Version:   1.4

##   Copyright (c) Richard Izzo, Luca Antiga, David Steinman. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

## Note: this class was contributed by
##       Richard Izzo (Github @rlizzo)
##       University at Buffalo
##       The Jacobs Institute

from __future__ import absolute_import #NEEDS TO STAY AS TOP LEVEL MODULE FOR Py2-3 COMPATIBILITY
import vtk
import sys

from vmtk import vtkvmtk
from . import vmtkrenderer
from vmtk import pypes
import pickle
import os

try:
    import numpy as np
except ImportError:
    raise ImportError('Unable to Import vmtknumpytosurface module, numpy is not installed')

class vmtkNumpyWriter(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.ArrayDict = None
        self.OutputFileName = ''
        self.Format = 'pickle'
        self.Compression = 1
        self.CompressionLevel = 4

        self.SetScriptName('vmtkNumpyWriter')
        self.SetScriptDoc('Writes a dictionary containing a nested dictionary of numpy arrays (generated by'
                          'vmtkcenterlinestonumpy, vmtkimagetonumpy, or vmtksurface to numpy) to disk as either'
                          'a python pickle object or as hdf5 file')

        self.SetInputMembers([
            ['ArrayDict','i','dict',1,'','the input array dictionary','vmtknumpyreader'],
            ['OutputFileName','ofile','str',1,'','the output file name'],
            ['Compression','compression','str',1,'(1,0)','Boolean value to compress hdf5 dataset files using gzip, default=1'],
            ['CompressionLevel','compressionlevel','str',1,'(0,9)','Specify compression level for gzip compressed hdf5 files. '
                                                                   'Must be an intiger from 0 to 9, higher levels have more compression,'
                                                                   'but take longer to process. Default=4'],
            ['Format','format','str',1,'["pickle","hdf5"]','write files as pickled object or hdf5 file format']])
        self.SetOutputMembers([])

    def WritePickledObjectFile(self):
        self.PrintLog('Writing Pickled Object File')
        pickleFileName = self.OutputFileName + '.pickle'
        with open(pickleFileName, 'wb') as outfile:
            pickle.dump(self.ArrayDict, outfile, protocol=pickle.HIGHEST_PROTOCOL)

        return

    def WriteHDF5File(self): # dic, filename, objectname=None):
        """
        Save a dictionary to an HDF5 file.
        """

        try:
            import h5py
        except ImportError:
            self.PrintError('ImportError: Unable to Write to hdf5. h5py module not installed')

        def recursively_save_dict_contents_to_group(h5file, path, dic):
            """
            Take an already open HDF5 file and insert the contents of a dictionary
            at the current path location. Can call itself recursively to fill
            out HDF5 files with the contents of a dictionary.
            """
            for key, item in dic.items():
                if isinstance(item, dict):
                    if not item.items():
                        h5file.create_group(path + key)
                    else:
                        recursively_save_dict_contents_to_group(h5file, path + key + '/', item)

                elif isinstance(item, list):
                    for index, element in enumerate(item):
                        if not isinstance(index, str):
                            index = str(index) 
                        if self.Compression == 1:
                            h5file.create_dataset(path + key + '/' + index, data=element, compression='gzip',
                                                  compression_opts=self.CompressionLevel)
                        else:
                            h5file[path + key + '/' + index] = element

                else:
                    if self.Compression == 1:
                        h5file.create_dataset(path+key, data=item, compression='gzip',
                                              compression_opts=self.CompressionLevel)
                    else:
                        h5file[path + key] = item

        self.PrintLog('Writing HDF5 File')
        hdf5FileName = self.OutputFileName + '.hdf5'
        with h5py.File(hdf5FileName, 'w') as h5file:
            recursively_save_dict_contents_to_group(h5file, '/', self.ArrayDict)

        return

    def Execute(self):

        if self.ArrayDict == None:
            self.PrintError('Error: no input dictionary')

        if self.OutputFileName == '':
            self.PrintError('Error: no output file name specified')

        if self.OutputFileName == 'BROWSER':
            import tkinter.filedialog
            import os.path
            initialDir = pypes.pypeScript.lastVisitedPath
            self.OutputFileName = tkinter.filedialog.asksaveasfilename(title="Output Dictionary",initialdir=initialDir)
            pypes.pypeScript.lastVisitedPath = os.path.dirname(self.OutputFileName)
            if not self.OutputFileName:
                self.PrintError('Error: no Output File Name.')

        self.OutputFileName = self.OutputFileName.rsplit( ".", 1 )[ 0 ]

        self.PrintLog('Writing File')
        if self.Format == 'pickle':
            self.WritePickledObjectFile()

        elif self.Format == 'hdf5':
            if self.Compression == 1:
                if not isinstance(self.CompressionLevel,int):
                    self.PrintError('Error, compression level must be of type Int')
                if self.CompressionLevel not in range(0, 10):
                    self.PrintError('Error, compression level '+ str(self.CompressionLevel) + ' is not in the valid range between'
                                                                                              '0 and 9 (inclusive).')
            self.WriteHDF5File()

        else:
            self.PrintError('Error: unsupported format '+ self.Format + '.')


if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()