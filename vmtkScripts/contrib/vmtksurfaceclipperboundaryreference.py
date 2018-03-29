#!/usr/bin/env python

## Program:   VMTK
## Module:    $RCSfile: surfaceclippercenterline2.py,v $
## Language:  Python
## Date:      $Date: 2017/10/09 11:34:08 $
## Version:   $Revision: 1.0 $

##   Copyright (c) Luca Antiga, David Steinman. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

## Note: this class was contributed by
##       Simone Manini
##       Orobix Srl
## And adapted and modified by
##       Kurt Sansom

from __future__ import absolute_import #NEEDS TO STAY AS TOP LEVEL MODULE FOR Py2-3 COMPATIBILITY
import vtk
import sys
from vmtk import vtkvmtk
from vmtk import vmtkrenderer
from vmtk import pypes
import copy


class vmtksubsetremesh():

    def __init__(self):

        self.Surface = None

        self.TargetArea = 1.0
        self.TargetEdgeLength = 1.0
        self.TargetAreaFactor = 1.0
        self.TargetEdgeLengthFactor = 1.0
        self.TriangleSplitFactor = 5.0
        self.MaxArea = 1E16
        self.MinArea = 0.0
        self.MaxEdgeLength = 1E16
        self.MinEdgeLength = 0.0
        self.NumberOfIterations = 10
        self.NumberOfConnectivityOptimizationIterations = 20
        self.CellEntityIdsArrayName = None
        self.TargetAreaArrayName = 'TargetArea'
        self.TargetEdgeLengthArrayName = ''
        self.ElementSizeMode = 'edgelength' #'area'
        self.MinAreaFactor = 0.5
        self.AspectRatioThreshold = 1.2
        self.InternalAngleTolerance = 0.0
        self.NormalAngleTolerance = 0.2
        self.CollapseAngleThreshold = 0.2
        self.Relaxation = 0.5
        self.PreserveBoundaryEdges = 0
        self.ExcludeEntityIds = []

    def Execute(self):

        if self.Surface == None:
            self.PrintError('Error: No input surface.')

        if self.ElementSizeMode == 'edgelength':
            self.TargetArea = 0.25 * 3.0**0.5 * self.TargetEdgeLength**2
        elif self.ElementSizeMode == 'edgelengtharray':
            calculator = vtk.vtkArrayCalculator()
            calculator.SetInputData(self.Surface)
            calculator.AddScalarArrayName(self.TargetEdgeLengthArrayName,0)
            calculator.SetFunction("%f^2 * 0.25 * sqrt(3) * %s^2" % (self.TargetEdgeLengthFactor,self.TargetEdgeLengthArrayName))
            calculator.SetResultArrayName(self.TargetAreaArrayName)
            calculator.Update()
            self.MaxArea = 0.25 * 3.0**0.5 * self.MaxEdgeLength**2
            self.MinArea = 0.25 * 3.0**0.5 * self.MinEdgeLength**2
            self.Surface = calculator.GetOutput()

        excludedIds = vtk.vtkIdList()
        if self.ExcludeEntityIds:
            for excludedId in self.ExcludeEntityIds:
                excludedIds.InsertNextId(excludedId)

        surfaceRemeshing = vtkvmtk.vtkvmtkPolyDataSurfaceRemeshing()
        surfaceRemeshing.SetInputData(self.Surface)
        if self.CellEntityIdsArrayName:
            surfaceRemeshing.SetCellEntityIdsArrayName(self.CellEntityIdsArrayName)
        if self.ElementSizeMode in ['area','edgelength']:
            surfaceRemeshing.SetElementSizeModeToTargetArea()
        elif self.ElementSizeMode in ['areaarray','edgelengtharray']:
            surfaceRemeshing.SetElementSizeModeToTargetAreaArray()
            surfaceRemeshing.SetTargetAreaArrayName(self.TargetAreaArrayName)
        else:
            self.PrintError('Error: unsupported ElementSizeMode.')
        surfaceRemeshing.SetTargetArea(self.TargetArea)
        surfaceRemeshing.SetTargetAreaFactor(self.TargetAreaFactor)
        surfaceRemeshing.SetTriangleSplitFactor(self.TriangleSplitFactor)
        surfaceRemeshing.SetMaxArea(self.MaxArea)
        surfaceRemeshing.SetMinArea(self.MinArea)
        surfaceRemeshing.SetNumberOfIterations(self.NumberOfIterations)
        surfaceRemeshing.SetNumberOfConnectivityOptimizationIterations(self.NumberOfConnectivityOptimizationIterations)
        surfaceRemeshing.SetRelaxation(self.Relaxation)
        surfaceRemeshing.SetMinAreaFactor(self.MinAreaFactor)
        surfaceRemeshing.SetAspectRatioThreshold(self.AspectRatioThreshold)
        surfaceRemeshing.SetInternalAngleTolerance(self.InternalAngleTolerance)
        surfaceRemeshing.SetNormalAngleTolerance(self.NormalAngleTolerance)
        surfaceRemeshing.SetCollapseAngleThreshold(self.CollapseAngleThreshold)
        surfaceRemeshing.SetPreserveBoundaryEdges(self.PreserveBoundaryEdges)
        surfaceRemeshing.SetExcludedEntityIds(excludedIds)
        surfaceRemeshing.Update()

        self.Surface = surfaceRemeshing.GetOutput()



class vmtkSurfaceClipperBoundaryReference(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.BoundaryInfo = None

        self.Surface = None

        self.CleanOutput = 1
        self.remesh = 1

        self.SetScriptName('vmtksurfaceclippercenterline')
        self.SetScriptDoc('interactively clip a tubular surface with normals estimated from centerline tangents')
        self.SetInputMembers([
            ['Surface','i','vtkPolyData',1,'','the input surface','vmtksurfacereader'],
            ['BoundaryInfo','boundaryinfo','vtkPolyData',1,'','the boundary info file','vmtksurfacereader'],
            ])
        self.SetOutputMembers([
            ['Surface','o','vtkPolyData',1,'','the output surface','vmtksurfacewriter'],
            ])

    def Execute(self):

        if self.Surface == None:
            self.PrintError('Error: no Surface.')

        if self.BoundaryInfo == None:
            self.PrintError('Error: No input Boundary reference system.')

        surfaceCleaner = vtk.vtkCleanPolyData()
        surfaceCleaner.SetInputData(self.Surface)
        surfaceCleaner.Update()

        surfaceTriangulator = vtk.vtkTriangleFilter()
        surfaceTriangulator.SetInputConnection(surfaceCleaner.GetOutputPort())
        surfaceTriangulator.PassLinesOff()
        surfaceTriangulator.PassVertsOff()
        surfaceTriangulator.Update()

        # Generate data arrays containing point and cell ids
        ids_name = "Ids"
        ids = vtk.vtkIdFilter()
        ids.SetInputConnection(surfaceTriangulator.GetOutputPort())
        #ids.PointIdsOn()
        ids.CellIdsOn()
        #ids.FieldDataOn()
        ids.SetIdsArrayName(ids_name)
        ids.Update()

        #writer2 = vtk.vtkXMLPolyDataWriter()
        #writer2.SetFileName("cell_ids_test.vtp")
        #writer2.SetInputConnection(ids.GetOutputPort())
        #writer2.Update()

        clippedSurface = ids.GetOutput()

       
        # Extract the boundary reference points from  polydata
        brd = self.BoundaryInfo
        # Get the number of points in the polydata
        idNumPointsInFile = brd.GetNumberOfPoints()
        narrays = brd.GetPointData().GetNumberOfArrays()
        for i in range(narrays):
            arrayName =  brd.GetPointData().GetArrayName(i)
            if(arrayName == "BoundaryNormals"):
                brd_normals = brd.GetPointData().GetArray(i)
        
        for i in range(idNumPointsInFile):
            br_point = brd.GetPoint(i)
            br_normal = brd_normals.GetTuple(i)
            inward_normal = tuple( -p for p in br_normal)

            locator = vtk.vtkPointLocator()
            locator.SetDataSet(clippedSurface)
            locator.BuildLocator()

            seedPointId = locator.FindClosestPoint(br_point)
            surf_point = clippedSurface.GetPoint(seedPointId)
            
            pt_diff = [0.0, 0.0, 0.0]
            vtk.vtkMath.Subtract(surf_point, br_point, pt_diff)
            pt_normalize = vtk.vtkMath.Normalize(pt_diff)
            
            ortho = vtk.vtkMath.Dot(inward_normal, pt_diff)
            #print(ortho)
            if (abs(ortho) > 0.49):
                #proceed
                continue
            
            plane = vtk.vtkPlane()
            plane.SetOrigin(br_point)
            plane.SetNormal(inward_normal)

            seamFilter = vtkvmtk.vtkvmtkTopologicalSeamFilter()
            seamFilter.SetInputData(clippedSurface)
            seamFilter.SetClosestPoint(surf_point)
            seamFilter.SetSeamScalarsArrayName("SeamScalars")
            seamFilter.SetSeamFunction(plane)
            seamFilter.Update()
            
            #print("got here")
            #writer3 = vtk.vtkXMLPolyDataWriter()
            #writer3.SetFileName("seam_test_{0}.vtp".format(i))
            #writer3.SetInputConnection(seamFilter.GetOutputPort())
            #writer3.Update()


            tree = vtk.vtkModifiedBSPTree()
            tree.SetDataSet(seamFilter.GetOutput())
            tree.BuildLocator()
            #intersect the locator with the line
            LineP0 = br_point
            dt = 0.1

            # 200 points
            dtheta = [0.0] #[0.0 + i*(359.0-0.0)/(300-1) for i in range(300)]
            out_vector = (0.0,0.0,0.0)
            tolerance = 0.0000001

            IntersectPointsList = vtk.vtkPoints()
            IntersectCellsList = vtk.vtkIdList()
            for theta in dtheta:
                IntersectPoints = vtk.vtkPoints()
                IntersectCells = vtk.vtkIdList()
                code = 0
                count = 1
                rotate = vtk.vtkTransform()
                rotate.RotateWXYZ(theta,  inward_normal)
                rotate.Update()

                #print(dir(rotate))
                #trans_m = vtk.vtkMatrix4x4()
                #rotate.GetMatrix(trans_m)

                out_vector = rotate.TransformVector(pt_diff)
                LineP1 = [ c2 + count*dt*c1 for c2, c1 in zip(br_point, out_vector)]
                #print(pt_diff, out_vector)
                while ( code == 0 and count < 10000):
                    count += 1
                    code = tree.IntersectWithLine(LineP0, LineP1,
                                                    tolerance, IntersectPoints,
                                                    IntersectCells)
                    LineP1 = [ c2 + count*dt*c1 for c2, c1 in zip(br_point, out_vector)]
                if(count > 10000 and code == 0):
                    print("no intersection")
                    continue

                if (code != 0):
                    pt = IntersectPoints.GetPoint(0)
                    #pt = [ c2 + dt*c1 for c2, c1 in zip(pt, out_vector)] # add some buffer, may not need it
                    IntersectPointsList.InsertNextPoint(pt)
                    IntersectCellsList.InsertNextId(IntersectCells.GetId(0))
                    #print(IntersectPoints.GetPoint(0), IntersectCells.GetId(0) )


                #Create a cell array to connect the points into meaningful geometry
                vertexIndices = vtk.vtkIdList()
                numPts = IntersectPointsList.GetNumberOfPoints()

                for i in range(numPts):
                    vertexIndices.InsertNextId(i)
                # Set the last vertex to 0; this means the last line segment will join the 19th point (vertices[19])
                # with the first one (vertices[0]), thus closing the circle.
                vertexIndices.InsertNextId(0)

                lines = vtk.vtkCellArray()
                lines.InsertNextCell(vertexIndices)
                # Create polydata to hold the geometry just created, and populate it
                test_data = vtk.vtkPolyData()
                test_data.SetPoints(IntersectPointsList)
                test_data.SetLines(lines);
                # writer4 = vtk.vtkXMLPolyDataWriter()
                # writer4.SetFileName("intersection_test_{0}.vtp".format(seed_count))
                # writer4.SetInputData(test_data)
                # writer4.Update()

                ### tried using threshold but to no avail
                thresh = vtk.vtkThreshold()
                thresh.SetInputConnection(seamFilter.GetOutputPort())
                thresh.ThresholdBetween(-0.99, 0.99)
                thresh.SetInputArrayToProcess(1, 0, 0, 0, "SeamScalars")
                thresh.Update()

                #pointIdArray = thresh.GetOutput().GetPointData().GetArray(ids_name)
                cellIdArray = thresh.GetOutput().GetCellData().GetArray(ids_name)
                #print(dir(cellIdArray))
                # Find all cells connected to point 0
                n_cells = cellIdArray.GetNumberOfTuples()
                print(n_cells)
                neighborhood = []
                visited_list = set()
                new_neighbors = set()
                cell_id_list = vtk.vtkIdTypeArray()
                cell_id_list.SetNumberOfComponents(1)

                neighbor_list = set()
                for c in range(n_cells):
                    neighbor_list.add(int(cellIdArray.GetTuple(c)[0]))

                # number of iterations to grow neighbors
                for i in range(3):
                    print("neighbor list", len(neighbor_list))
                    for cell_id in neighbor_list:
                        #print(cell_id)
                        visited_list.add(cell_id)
                        cellPointIds = vtk.vtkIdList()
                        seamFilter.GetOutput().GetCellPoints(cell_id, cellPointIds)
                        #neighbor cells may be listed multiple times
                        # set instead of list if you want a unique list of neighbors

                        # For each vertice of the cell, we calculate which cells uses that point.
                        # So if we make this, for each vertice, we have all the neighbors.
                        # In the case we use ''cellPointIds'' as a parameter of ''GeteCellNeighbors'',
                        # we will obtain an empty set. Because the only cell that is using that set of points
                        # is the current one. That is why we have to make each vertice at time.
                        for i in range(cellPointIds.GetNumberOfIds()):
                            idList = vtk.vtkIdList()
                            idList.InsertNextId(cellPointIds.GetId(i))
                            #get the neighbors of the cell
                            neighborCellIds = vtk.vtkIdList()
                            seamFilter.GetOutput().GetCellNeighbors(cell_id, idList, neighborCellIds)
                            for j in range(neighborCellIds.GetNumberOfIds()):
                                new_neighbors.add(neighborCellIds.GetId(j))
                    # get unique neighbors
                    neighbor_list = new_neighbors.difference(visited_list)
                    new_neighbors.clear()

                # unique list of cell Ids
                hood_list = list(visited_list.union(neighbor_list)) # list(visited_list.union(new_neighbors))

                ### tried remeshing subset but doesn't connect with parent
                for idx in hood_list:
                    cell_id_list.InsertNextValue(idx)

                cell_entity_id_array = vtk.vtkUnsignedCharArray()
                cell_entity_id_array.SetName("RemeshCells")
                cell_entity_id_array.SetNumberOfTuples(seamFilter.GetOutput().GetNumberOfCells())
                cell_entity_id_array.FillComponent(0,0)

                for cell_idx in list(visited_list):
                    cell_entity_id_array.SetValue(cell_idx, 1)

                seamFilter.GetOutput().GetCellData().AddArray(cell_entity_id_array)
                seamFilter.Update()

                selectionNode = vtk.vtkSelectionNode()
                selectionNode.SetFieldType(vtk.vtkSelectionNode().CELL)
                selectionNode.SetContentType(vtk.vtkSelectionNode().INDICES)
                selectionNode.SetSelectionList(cell_id_list)

                selection = vtk.vtkSelection()
                selection.AddNode(selectionNode)

                extractSelectedIds = vtk.vtkExtractSelectedIds()
                extractSelectedIds.SetInputData(0, seamFilter.GetOutput())
                extractSelectedIds.SetInputData(1, selection);
                extractSelectedIds.Update()

                geometryFilter = vtk.vtkDataSetSurfaceFilter() #vtkGeometryFilter()
                geometryFilter.SetInputConnection(extractSelectedIds.GetOutputPort())
                geometryFilter.Update()
                geometryFilter.GetOutput().BuildLinks()

                dist_2_plane = vtk.vtkDoubleArray()
                dist_2_plane.SetName("PlaneDistance")
                dist_2_plane.SetNumberOfComponents(1)
                dist_2_plane.SetNumberOfTuples(geometryFilter.GetOutput().GetNumberOfPoints())
                dist_2_plane.FillComponent(0, 0.0)

                for pt_id in range(geometryFilter.GetOutput().GetNumberOfPoints()):
                    pt = geometryFilter.GetOutput().GetPoint(pt_id)
                    dist = abs(plane.FunctionValue(pt))/2.0
                    if (dist > 0.5):
                        dist = 0.5
                    elif (dist < 0.02):
                        dist = 0.02
                    dist_2_plane.SetValue(pt_id, dist)

                geometryFilter.GetOutput().GetPointData().AddArray(dist_2_plane)
                geometryFilter.Update()


                remeshsubset = vmtksubsetremesh()
                remeshsubset.CellEntityIdsArrayName = "RemeshCells"
                remeshsubset.ElementSizeMode = 'edgelengtharray'
                remeshsubset.TargetEdgeLengthArrayName = "PlaneDistance"
                remeshsubset.TargetEdgeLengthFactor = 0.8
                remeshsubset.NumberOfIterations = 10
                remeshsubset.Surface = geometryFilter.GetOutput()
                remeshsubset.ExcludeEntityIds = [0] #list(free_edge_cells)
                remeshsubset.Execute()


                #invert the selection
                selectionNode.GetProperties().Set(vtk.vtkSelectionNode().INVERSE(), 1)
                extractSelectedIds.Update()
                geometryFilter.Update()

                #Append the two meshes
                appendFilter = vtk.vtkAppendPolyData()
                appendFilter.AddInputData(remeshsubset.Surface)
                appendFilter.AddInputData(geometryFilter.GetOutput())
                appendFilter.Update()

                # Remove any duplicate points.
                cleanFilter = vtk.vtkCleanPolyData()
                cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
                cleanFilter.Update()

                surfaceTriangulator2 = vtk.vtkTriangleFilter()
                surfaceTriangulator2.SetInputConnection(cleanFilter.GetOutputPort())
                surfaceTriangulator2.PassLinesOff()
                surfaceTriangulator2.PassVertsOff()
                surfaceTriangulator2.Update()


                seamFilter2 = vtkvmtk.vtkvmtkTopologicalSeamFilter()
                seamFilter2.SetInputConnection(surfaceTriangulator2.GetOutputPort())
                seamFilter2.SetClosestPoint(IntersectPointsList.GetPoint(0))
                seamFilter2.SetSeamScalarsArrayName("SeamScalars")
                seamFilter2.SetSeamFunction(plane)
                seamFilter2.Update()


                clipper = vtk.vtkClipPolyData()
                clipper.SetInputConnection(seamFilter2.GetOutputPort())
                #clipper.SetClipFunction(plane)
                clipper.GenerateClipScalarsOff()
                clipper.GenerateClippedOutputOn()
                #clipper.InsideOutOff()
                clipper.Update()

                #writer7 = vtk.vtkXMLPolyDataWriter()
                #writer7.SetFileName("clipper_test{0}.vtp".format(i))
                #writer7.SetInputConnection(clipper.GetOutputPort())
                #writer7.Update()

                connectivity = vtk.vtkPolyDataConnectivityFilter()
                connectivity.SetInputConnection(clipper.GetOutputPort())
                connectivity.SetExtractionModeToLargestRegion()
                #connectivity.SetExtractionModeToClosestPointRegion()
                #connectivity.SetClosestPoint(mainBodyPoint)
                connectivity.Update()

                surfaceCleaner = vtk.vtkCleanPolyData()
                surfaceCleaner.SetInputConnection(connectivity.GetOutputPort())
                surfaceCleaner.Update()

                surfaceTriangulator = vtk.vtkTriangleFilter()
                surfaceTriangulator.SetInputConnection(surfaceCleaner.GetOutputPort())
                surfaceTriangulator.PassLinesOff()
                surfaceTriangulator.PassVertsOff()
                surfaceTriangulator.Update()

                ids_filter = vtk.vtkIdFilter()
                ids_filter.SetInputConnection(surfaceTriangulator.GetOutputPort())
                #ids_filter.PointIdsOn()
                ids_filter.CellIdsOn()
                #ids_filter.FieldDataOn()
                ids_filter.SetIdsArrayName(ids_name)
                ids_filter.Update()

                clippedSurface = vtk.vtkPolyData()
                clippedSurface.DeepCopy(ids_filter.GetOutput())
                #clippedSurface = surfaceTriangulator.GetOutput()

        self.Surface = clippedSurface


if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
