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


class vmtkFlowExtensions2(pypes.pypeScript):

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
        n_clip_pts = reader.GetOutput().GetNumberOfPoints()
        print("nyumber of points: ", n_clip_pts)

        locator_ctrpt = vtk.vtkPointLocator()
        locator_ctrpt.SetDataSet(seedPolyData)
        locator_ctrpt.BuildLocator()
        case_dict = {"inlet" : vtk.vtkIdList(), "outlet" : vtk.vtkIdList()}
        for clip_pt_id in range(n_clip_pts):
            clip_pt = reader.GetOutput().GetPoint(clip_pt_id)
            # 0 for inlet and 1 for outlet
            flow_dir = reader.GetOutput().GetCellData().GetArray("BoundaryRegion").GetTuple(clip_pt_id)

            bary_pt_id = locator_ctrpt.FindClosestPoint(clip_pt)
            bary_pt = seedPolyData.GetPoint(bary_pt_id)

            dist = (vtk.vtkMath.Distance2BetweenPoints(clip_pt, bary_pt))**0.5
            #print("distance between points: ", dist)
            if (dist > 5.0):
                continue
            if (flow_dir[0] < 1.0):
                case_dict["inlet"].InsertNextId(bary_pt_id)
            else:
                case_dict["outlet"].InsertNextId(bary_pt_id)

        ctr_points = vtk.vtkPoints()
        ctr_points_id = vtk.vtkUnsignedCharArray()
        ctr_points_id.SetNumberOfComponents(1)
        ctr_points_id.SetName("BoundaryRegion")
        vertices = vtk.vtkCellArray()
        new_vertex_id = 0
        source_list = []
        target_list = []

        for (d, bndIds) in case_dict.items():
            flowExtensionsFilter = vtkvmtk.vtkvmtkPolyDataFlowExtensionsFilter()
            flowExtensionsFilter.SetInputData(extendedsurface)
            flowExtensionsFilter.SetCenterlines(self.Centerlines)
            flowExtensionsFilter.SetSigma(self.Sigma)
            flowExtensionsFilter.SetAdaptiveExtensionLength(self.AdaptiveExtensionLength)
            flowExtensionsFilter.SetAdaptiveExtensionRadius(self.AdaptiveExtensionRadius)
            flowExtensionsFilter.SetAdaptiveNumberOfBoundaryPoints(self.AdaptiveNumberOfBoundaryPoints)
            flowExtensionsFilter.SetExtensionLength(self.ExtensionLength)
            flowExtensionsFilter.SetExtensionRadius(self.ExtensionRadius)
            flowExtensionsFilter.SetCenterlineNormalEstimationDistanceRatio(self.CenterlineNormalEstimationDistanceRatio)
            flowExtensionsFilter.SetNumberOfBoundaryPoints(self.TargetNumberOfBoundaryPoints)
            flowExtensionsFilter.SetExtensionModeToUseCenterlineDirection()
            flowExtensionsFilter.SetBoundaryIds(bndIds)

            if( d == "inlet" ): #inlet
                flowExtensionsFilter.SetExtensionRatio(4.0)
                flowExtensionsFilter.SetTransitionRatio(0.25)
            else:
                flowExtensionsFilter.SetExtensionRatio(2.0)
                flowExtensionsFilter.SetTransitionRatio(0.75)

            flowExtensionsFilter.Update()
            extendedsurface = flowExtensionsFilter.GetOutput()

            # keep information about the inlets and outlets
            new_boundaryIds = vtk.vtkIdList()

            bndryExtractor = vtkvmtk.vtkvmtkPolyDataBoundaryExtractor()
            bndryExtractor.SetInputData(extendedsurface)
            bndryExtractor.Update()
            bndries = bndryExtractor.GetOutput()
            n_bndries = bndries.GetNumberOfCells()
            new_seedPoints = vtk.vtkPoints()
            for i in range(n_bndries):
                new_barycenter = [0.0, 0.0, 0.0]
                vtkvmtk.vtkvmtkBoundaryReferenceSystems.ComputeBoundaryBarycenter(bndries.GetCell(i).GetPoints(),new_barycenter)

                bary_pt_id = locator_ctrpt.FindClosestPoint(new_barycenter)
                bary_pt = seedPolyData.GetPoint(bary_pt_id)

                dist = (vtk.vtkMath.Distance2BetweenPoints(new_barycenter, bary_pt))**0.5
                #print("distance between points: ", dist)
                if (dist < 1.0E-4):
                    continue
                else:
                    ctr_points.InsertNextPoint(new_barycenter)
                    vertex = vtk.vtkVertex()
                    vertex.GetPointIds().InsertNextId(int(new_vertex_id))
                    new_vertex_id += 1
                    vertices.InsertNextCell(vertex)

                    if( d == "outlet"):
                        ctr_points_id.InsertNextTuple([int(1)])
                        target_list.append(new_barycenter)
                    else:
                        ctr_points_id.InsertNextTuple([int(0)])
                        source_list.append(new_barycenter)

                    seedPolyData.GetPoints().SetPoint(bary_pt_id, new_barycenter)
                    locator_ctrpt.SetDataSet(seedPolyData)
                    locator_ctrpt.BuildLocator()



        surfaceCapper = vtkvmtk.vtkvmtkCapPolyData()
        surfaceCapper.SetInputData(extendedsurface)
        surfaceCapper.SetDisplacement(0.0)
        surfaceCapper.SetInPlaneDisplacement(0.0)
        surfaceCapper.Update()
        centerlineInputSurface = surfaceCapper.GetOutput()
        capCenterIds = surfaceCapper.GetCapCenterIds()
        capBoundaryIds = surfaceCapper.GetBoundaryIds()

        boundary_array = vtk.vtkUnsignedCharArray()
        boundary_array.SetNumberOfComponents(1)
        boundary_array.SetNumberOfTuples(capCenterIds.GetNumberOfIds())
        boundary_array.SetName("BoundaryId")
        #print(capCenterIds, dir(capCenterIds))
        print(capCenterIds.GetNumberOfIds())

        seedPolyData.SetPoints(ctr_points)
        locator_ctrpt.SetDataSet(seedPolyData)
        locator_ctrpt.BuildLocator()
        #print(dir(capBoundaryIds))
        for i in range(capCenterIds.GetNumberOfIds()):
            ctr_pt_id = capCenterIds.GetId(i)
            center_pt = centerlineInputSurface.GetPoint(ctr_pt_id)
            bary_pt_id = locator_ctrpt.FindClosestPoint(center_pt)
            bary_pt = seedPolyData.GetPoint(bary_pt_id)

            dist = (vtk.vtkMath.Distance2BetweenPoints(center_pt, bary_pt))**0.5
            print(" distance between: ", dist)

            #bndId = capBoundaryIds.GetId(ctr_pt_id)
            flow_direction = ctr_points_id.GetTuple(bary_pt_id)
            print("stuff: ", i, flow_direction[0])
            boundary_array.SetTuple(bary_pt_id, [int(i+2)])

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(ctr_points)
        polydata.SetVerts(vertices)
        polydata.GetCellData().AddArray(ctr_points_id)
        polydata.GetCellData().AddArray(boundary_array)
        print(ctr_points.GetNumberOfPoints(), ctr_points_id.GetNumberOfTuples())

        writer7 = vtk.vtkXMLPolyDataWriter()
        writer7.SetFileName("centerline_points_ext.vtp")
        writer7.SetInputData(polydata)
        writer7.Update()

        #print("source pts: ", [item for sublist in source_list for item in sublist])
        #print("target pts: ", [item for sublist in target_list for item in sublist])

        print("source pts: ", " ".join(map(str, [item for sublist in source_list for item in sublist])))
        print("target pts: ", " ".join(map(str, [item for sublist in target_list for item in sublist])))
        self.Surface = extendedsurface



if __name__=='__main__':

    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
