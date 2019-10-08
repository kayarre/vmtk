#!/usr/bin/env python

## Program:   VMTK
## Module:    $RCSfile: vmtkmesharrayoperation.py,v $
## Language:  Python
## Date:      $Date: 2005/09/14 09:49:59 $
## Version:   $Revision: 1.7 $

##   Copyright (c) Luca Antiga, David Steinman. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

from __future__ import absolute_import #NEEDS TO STAY AS TOP LEVEL MODULE FOR Py2-3 COMPATIBILITY
import vtk
from vmtk import vtkvmtk
import sys

from vmtk import pypes


class vmtkMeshArrayOperation(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)
        
        self.Mesh = None
        self.Mesh2 = None
        self.InputArrayName = None
        self.Input2ArrayName = None
        self.ResultArrayName = "Result"
        self.Constant = 0.0

        self.Operation = 'add'

        self.SetScriptName('vmtkmesharrayoperation')
        self.SetScriptDoc('perform an operation between arrays of two meshes and store result in the first mesh')
        self.SetInputMembers([
            ['Mesh','i','vtkUnstructuredGrid',1,'','the input mesh','vmtkmeshreader'],
            ['Mesh2','i2','vtkUnstructuredGrid',1,'','the second input mesh','vmtkmeshreader'],
            ['Operation','operation','str',1,'["multiplybyc","addc","add","subtract","min","max"]','the operation to be performed on the array; multiplybyc and addc only require the first input mesh to be specified'],
            ['Constant','constant','float',1,'','the value of the constant for multiplybyc and addc'],
            ['InputArrayName','iarray','str',1,'','the name of the array on the first mesh'],
            ['Input2ArrayName','i2array','str',1,'','the name of the array on the second mesh; if unspecified, InputArrayName is used'],
            ['ResultArrayName','resultarray','str',1,'','the name of the array where the result of the operation is stored']
            ])
        self.SetOutputMembers([
            ['Mesh','o','vtkUnstructuredGrid',1,'','the output mesh','vmtkmeshwriter']
            ])

    def Execute(self):

        if self.Mesh == None:
            self.PrintError('Error: No Mesh.')

        if self.Mesh2 == None and self.Operation not in ["multiplybyc", "addc"]:
            self.PrintError('Error: No Mesh2.')

        array1 = self.Mesh.GetPointData().GetArray(self.InputArrayName)

        array2 = None

        if self.Operation not in ["multiplybyc", "addc"]:
            if self.Input2ArrayName:
                array2 = self.Mesh2.GetPointData().GetArray(self.Input2ArrayName)
            else:
                array2 = self.Mesh2.GetPointData().GetArray(self.InputArrayName)

        resultArray = vtk.vtkDoubleArray()
        resultArray.DeepCopy(array1)
        resultArray.SetName(self.ResultArrayName)

        for i in range(array1.GetNumberOfTuples()):
            value1 = array1.GetTuple1(i)
            value2 = None
            if array2:
                value2 = array2.GetTuple1(i)
            resultValue = 0.0
            if self.Operation == "multiplybyc":
                resultValue = value1 * self.Constant
            elif self.Operation == "addc":
                resultValue = value1 + self.Constant
            elif self.Operation == "add":
                resultValue = value1 + value2
            elif self.Operation == "subtract":
                resultValue = value1 - value2
            elif self.Operation == "min":
                resultValue = min(value1, value2)
            elif self.Operation == "max":
                resultValue = max(value1, value2)
            resultArray.SetValue(i,resultValue)

        self.Mesh.GetPointData().AddArray(resultArray)


if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()

