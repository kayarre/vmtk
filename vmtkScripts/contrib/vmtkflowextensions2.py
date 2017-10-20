#!/usr/bin/env python

## Program:   VMTK
## Module:    $RCSfile: vmtkflowextensions2.py,v $
## Language:  Python
## Date:      $Date: 2006/07/17 09:52:56 $
## Version:   $Revision: 1.7 $

##   Copyright (c) Luca Antiga, David Steinman. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

from __future__ import absolute_import #NEEDS TO STAY AS TOP LEVEL MODULE FOR Py2-3 COMPATIBILITY
import vtk
import sys

from vmtk import vtkvmtk
from vmtk import pypes
from vmtk import vmtkrenderer


class vmtkFlowExtensions(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.Surface = None
        self.Centerlines = None

        self.AdaptiveExtensionLength = 0
        self.AdaptiveExtensionRadius = 1
        self.AdaptiveNumberOfBoundaryPoints = 0
        self.ExtensionLength = 1.0
        self.ExtensionRatio = 10.0
        self.ExtensionRadius = 1.0
        self.TransitionRatio = 0.25
        self.TargetNumberOfBoundaryPoints = 50
        self.CenterlineNormalEstimationDistanceRatio = 1.0
        self.ExtensionMode = "centerlinedirection"
        self.InterpolationMode = "thinplatespline"
        self.Sigma = 1.0

        self.SetScriptName('vmtkflowextensions')
        self.SetInputMembers([
            ['Surface','i','vtkPolyData',1,'','','vmtksurfacereader'],
            ['Centerlines','centerlines','vtkPolyData',1,'','','vmtksurfacereader'],
            ['ExtensionMode','extensionmode','str',1,'["centerlinedirection","boundarynormal"]','method for computing the normal for extension'],
            ['InterpolationMode','interpolationmode','str',1,'["linear","thinplatespline"]','method for computing interpolation from the model section to a circular section'],
            ['Sigma','sigma','float',1,'(0.0,)','thin plate spline stiffness'],
            ['AdaptiveExtensionLength','adaptivelength','bool',1],
            ['AdaptiveExtensionRadius','adaptiveradius','bool',1],
            ['AdaptiveNumberOfBoundaryPoints','adaptivepoints','bool',1],
            ['ExtensionLength','extensionlength','float',1,'(0.0,)'],
            ['ExtensionRatio','extensionratio','float',1,'(0.0,)'],
            ['ExtensionRadius','extensionradius','float',1,'(0.0,)'],
            ['TransitionRatio','transitionratio','float',1,'(0.0,)'],
            ['TargetNumberOfBoundaryPoints','boundarypoints','int',1,'(0,)'],
            ['CenterlineNormalEstimationDistanceRatio','normalestimationratio','float',1,'(0.0,)']
            ])
        self.SetOutputMembers([
            ['Surface','o','vtkPolyData',1,'','','vmtksurfacewriter'],
            ['Centerlines','centerlines','vtkPolyData',1]
            ])

    def LabelValidator(self,text):
        import string
        if not text:
            return 0
        if not text.split():
            return 0
        for char in text:
            if char not in string.digits + " ":
                return 0
        return 1

    def Execute(self):

        if self.Surface == None:
            self.PrintError('Error: No input surface.')

        if self.ExtensionMode == "centerlinedirection" and self.Centerlines == None:
            self.PrintError('Error: No input centerlines.')

        extendedsurface = vtk.vtkPolyData()
        extendedsurface.DeepCopy(self.Surface)

        boundaryIds = vtk.vtkIdList()

        boundaryExtractor = vtkvmtk.vtkvmtkPolyDataBoundaryExtractor()
        boundaryExtractor.SetInputData(self.Surface)
        boundaryExtractor.Update()
        boundaries = boundaryExtractor.GetOutput()
        numberOfBoundaries = boundaries.GetNumberOfCells()
        seedPoints = vtk.vtkPoints()
        for i in range(numberOfBoundaries):
            barycenter = [0.0, 0.0, 0.0]
            vtkvmtk.vtkvmtkBoundaryReferenceSystems.ComputeBoundaryBarycenter(boundaries.GetCell(i).GetPoints(),barycenter)
            seedPoints.InsertNextPoint(barycenter)
            boundaryIds.InsertNextId(i)
        seedPolyData = vtk.vtkPolyData()
        seedPolyData.SetPoints(seedPoints)

        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName("centerline_points.vtp")
        reader.Update()
        clip_pts = reader.GetOutput().GetNumberOfPoints()

        locator_ctrpt = vtk.vtkPointLocator()
        locator_ctrpt.SetDataSet(seedPolyData)
        locator_ctrpt.BuildLocator()
        case_dict = {"inlet" : vtk.vtkIdList(), "outlet" : vtk.vtkIdList()}
        for clip_pt_id in range(clip_pts):
            clip_pt = reader.GetOutput().GetPoint(clip_pt_id)
            # 0 for inlet and 1 for outlet
            flow_dir = reader.GetOutput().GetArray("BoundaryRegion").GetTuple(clip_pt)

            bary_pt_id = locator_ctrpt.FindClosestPoint(clip_pt)
            bary_pt = seedPolyData.GetPoint(bary_pt_id)

            dist = (vtk.vtkMath.Distance2BetweenPoints(clip_pt, bary_pt))**0.5
            print("distance between points ")
            if (dist > 5.0):
                continue
            if (flow_dir == 0):
                case_dict["inlet"].InsertNextId(bary_pt_id)
            else:
                case_dict["outlet"].InsertNextId(bary_pt_id)

        for (d, ob) in case_dict.items():
            flowExtensionsFilter = vtkvmtk.vtkvmtkPolyDataFlowExtensionsFilter()
            flowExtensionsFilter.SetInputData(extendedsurface)
            flowExtensionsFilter.SetCenterlines(self.Centerlines)
            flowExtensionsFilter.SetSigma(1.0)
            flowExtensionsFilter.SetAdaptiveExtensionLength(1)
            flowExtensionsFilter.SetAdaptiveExtensionRadius(0)
            flowExtensionsFilter.SetAdaptiveNumberOfBoundaryPoints(0)
            flowExtensionsFilter.SetExtensionLength(1.0)
            flowExtensionsFilter.SetExtensionRadius(1.0)
            flowExtensionsFilter.SetCenterlineNormalEstimationDistanceRatio(1.0)
            flowExtensionsFilter.SetNumberOfBoundaryPoints(200)
            flowExtensionsFilter.SetExtensionModeToUseCenterlineDirection()
            flowExtensionsFilter.SetBoundaryIds(boundaryIds)

            if( d == 0 ): #inlet
                flowExtensionsFilter.SetExtensionRatio(4.0)
                flowExtensionsFilter.SetTransitionRatio(0.25)
            else:
                flowExtensionsFilter.SetExtensionRatio(2.0)
                flowExtensionsFilter.SetTransitionRatio(0.75)

            flowExtensionsFilter.Update()
            extendedsurface = flowExtensionsFilter.GetOutput()

        self.Surface = extendedsurface



if __name__=='__main__':

    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
