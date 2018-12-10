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



## TODO: make SeedSelector a separate pype script to be used in other contexts
class vmtkSeedSelector(object):

    def __init__(self):
        self._Surface = None
        self._SeedIds = None
        self._MainBodySeedIds = vtk.vtkIdList()
        self._TargetSeedIds = vtk.vtkIdList()
        self._OutletSeedIds = vtk.vtkIdList()
        self.PrintError = None
        self.PrintLog = None
        self.InputText = None
        self.OutputText = None
        self.InputInfo = None

    def SetSurface(self,surface):
        self._Surface = surface
    def GetSurface(self):
        return self._Surface

    def GetMainBodySeedIds(self):
        return self._MainBodySeedIds

    def GetTargetSeedIds(self):
        return self._TargetSeedIds

    def GetOutletSeedIds(self):
        return self._OutletSeedIds

    def Execute(self):
        pass

class vmtkPickPointSeedSelector(vmtkSeedSelector):

    def __init__(self):
        vmtkSeedSelector.__init__(self)
        self.PickedSeedIds = vtk.vtkIdList()
        self.PickedSeeds = vtk.vtkPolyData()
        self.vmtkRenderer = None
        self.OwnRenderer = 0
        self.Script = None
        self.Opacity = 1
        self.continue_ = 0

    def UndoCallback(self, obj):
        self.InitializeSeeds()
        self.PickedSeeds.Modified()
        self.vmtkRenderer.RenderWindow.Render()
        
    def ContinueCallback(self, obj):
        self.continue_ = 1

    def PickCallback(self, obj):
        picker = vtk.vtkCellPicker()
        picker.SetTolerance(1E-4 * self._Surface.GetLength())
        eventPosition = self.vmtkRenderer.RenderWindowInteractor.GetEventPosition()
        result = picker.Pick(float(eventPosition[0]),float(eventPosition[1]),0.0,self.vmtkRenderer.Renderer)
        if result == 0:
            return
        pickPosition = picker.GetPickPosition()
        pickedCellPointIds = self._Surface.GetCell(picker.GetCellId()).GetPointIds()
        minDistance = 1E10
        pickedSeedId = -1
        for i in range(pickedCellPointIds.GetNumberOfIds()):
            distance = vtk.vtkMath.Distance2BetweenPoints(pickPosition,self._Surface.GetPoint(pickedCellPointIds.GetId(i)))
            if distance < minDistance:
                minDistance = distance
                pickedSeedId = pickedCellPointIds.GetId(i)
        if pickedSeedId == -1:
            pickedSeedId = pickedCellPointIds.GetId(0)
        self.PickedSeedIds.InsertNextId(pickedSeedId)
        point = self._Surface.GetPoint(pickedSeedId)
        self.PickedSeeds.GetPoints().InsertNextPoint(point)
        self.PickedSeeds.Modified()
        self.vmtkRenderer.RenderWindow.Render()

    def InitializeSeeds(self):
        self.PickedSeedIds.Initialize()
        self.PickedSeeds.Initialize()
        seedPoints = vtk.vtkPoints()
        self.PickedSeeds.SetPoints(seedPoints)

    def Execute(self):

        if (self._Surface == None):
            self.PrintError('vmtkPickPointSeedSelector Error: Surface not set.')
            return

        self._MainBodySeedIds.Initialize()
        self._TargetSeedIds.Initialize()
        self._OutletSeedIds.Initialize()

        if not self.vmtkRenderer:
            self.vmtkRenderer = vmtkrenderer.vmtkRenderer()
            self.vmtkRenderer.Initialize()
            self.OwnRenderer = 1

        self.vmtkRenderer.RegisterScript(self.Script)

        glyphs = vtk.vtkGlyph3D()
        glyphSource = vtk.vtkSphereSource()
        glyphs.SetInputData(self.PickedSeeds)
        glyphs.SetSourceConnection(glyphSource.GetOutputPort())
        glyphs.SetScaleModeToDataScalingOff()
        glyphs.SetScaleFactor(self._Surface.GetLength()*0.01)
        glyphMapper = vtk.vtkPolyDataMapper()
        glyphMapper.SetInputConnection(glyphs.GetOutputPort())
        self.SeedActor = vtk.vtkActor()
        self.SeedActor.SetMapper(glyphMapper)
        self.SeedActor.GetProperty().SetColor(1.0,0.0,0.0)
        self.SeedActor.PickableOff()
        self.vmtkRenderer.Renderer.AddActor(self.SeedActor)

        ##self.vmtkRenderer.RenderWindowInteractor.AddObserver("KeyPressEvent", self.KeyPressed)
        self.vmtkRenderer.AddKeyBinding('u','Undo.',self.UndoCallback)
        self.vmtkRenderer.AddKeyBinding('c','continue.',self.ContinueCallback)
        self.vmtkRenderer.AddKeyBinding('space','Add points.',self.PickCallback)

        surfaceMapper = vtk.vtkPolyDataMapper()
        surfaceMapper.SetInputData(self._Surface)
        surfaceMapper.ScalarVisibilityOff()
        surfaceActor = vtk.vtkActor()
        surfaceActor.SetMapper(surfaceMapper)
        surfaceActor.GetProperty().SetOpacity(1.0)

        self.vmtkRenderer.Renderer.AddActor(surfaceActor)

        self.InputInfo('Please position the mouse and press space to add a point on the main body of the model, \'u\' to undo\n')

        any = 0
        while any == 0:
            self.InitializeSeeds()
            self.vmtkRenderer.Render()
            any = self.PickedSeedIds.GetNumberOfIds()
        self._MainBodySeedIds.DeepCopy(self.PickedSeedIds)

        self.InputInfo('Please position the mouse and press space to add point on vessel wall near each Inlet location, \'u\' to undo\n')

        any = 0
        continue_ = 0
        while (any == 0) and (continue_ == 0):
            self.InitializeSeeds()
            self.vmtkRenderer.Render()
            any = self.PickedSeedIds.GetNumberOfIds()
            continue_ = self.continue_
        self._TargetSeedIds.DeepCopy(self.PickedSeedIds)

        self.InputInfo('Please position the mouse and press space to add point on vessel wall near each Outlet location, \'u\' to undo\n')

        any = 0
        continue_ = 0
        while (any == 0) and (continue_ == 0):
            self.InitializeSeeds()
            self.vmtkRenderer.Render()
            any = self.PickedSeedIds.GetNumberOfIds()
            continue_ = self.continue_
        self._OutletSeedIds.DeepCopy(self.PickedSeedIds)

        if self.OwnRenderer:
            self.vmtkRenderer.Deallocate()

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



class vmtkSurfaceClipperCenterline2(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.Centerlines = None
        self.BoundaryReference = None
        self.FrenetTangentArrayName = 'FrenetTangent'
        self.FrenetNormalArrayName = 'FrenetNormal'

        self.Surface = None
        self.vmtkRenderer = None
        self.OwnRenderer = 0
        self.SeedSelector = None
        self.SeedSelectorName = 'pointlist'

        self.SourcePoints = []
        self.TargetPoints = []
        self.MainBodyPointId = None
        
        self.Interactive = 1

        self.CleanOutput = 1
        self.remesh = 1


        self.SetScriptName('vmtksurfaceclippercenterline2')
        self.SetScriptDoc('interactively clip a tubular surface with normals estimated from centerline tangents')
        self.SetInputMembers([
            ['Surface','i','vtkPolyData',1,'','the input surface','vmtksurfacereader'],
            ['Centerlines','centerlines','vtkPolyData',1,'','the input centerlines','vmtksurfacereader'],
            ['BoundaryReference','boundaryreference','vtkPolyData',1,'','the input boundary reference file', 'vmtksurfacereader'],
            ['FrenetTangentArrayName','frenettangentarray','str',1,'','name of the array where centerline tangent vectors of the Frenet reference system are stored'],
            ['FrenetNormalArrayName','frenetnormalarray','str',1,'','name of the array where centerline normal vectors of the Frenet reference system are stored'],
            ['Interactive','interactive','bool',1],
            ['MainBodyPointId', 'mainbodypointid', 'int',1,'(0.0,)'],
            ['vmtkRenderer','renderer','vmtkRenderer',1,'','external renderer']
            ])
        self.SetOutputMembers([
            ['Surface','o','vtkPolyData',1,'','the output surface','vmtksurfacewriter'],
             ['BoundaryDirection','boundary','vtkPolyData',1,'','the vertex inlet/outlet point information','vmtksurfacewriter']
            ])

    def Execute(self):

        if self.Surface == None:
            self.PrintError('Error: no Surface.')

        if self.Centerlines == None:
            self.PrintError('Error: No input centerlines.')

        if self.Interactive:
            
            if not self.vmtkRenderer:
                self.vmtkRenderer = vmtkrenderer.vmtkRenderer()
                self.vmtkRenderer.Initialize()
                self.OwnRenderer = 1

            self.SeedSelector = vmtkPickPointSeedSelector()
            self.SeedSelector.vmtkRenderer = self.vmtkRenderer
            self.SeedSelector.Script = self

            self.SeedSelector.SetSurface(self.Surface)
            self.SeedSelector.InputInfo = self.InputInfo
            self.SeedSelector.InputText = self.InputText
            self.SeedSelector.OutputText = self.OutputText
            self.SeedSelector.PrintError = self.PrintError
            self.SeedSelector.PrintLog = self.PrintLog
            self.SeedSelector.Execute()

            mainBodySeedIds = self.SeedSelector.GetMainBodySeedIds()
            mainBodySeedId = mainBodySeedIds.GetId(0)

            targetSeedIds = self.SeedSelector.GetTargetSeedIds()
            outletSeedIds = self.SeedSelector.GetOutletSeedIds()
        
        else:
            if self.BoundaryReference == None:
                self.PrintError('Error: No Boundary reference system to define clipping location.')
            if self.MainBodyPointId == None:
                self.PrintError('Error: No main body point specified with -mainbodypointid.')
            
            locator_pt = vtk.vtkPointLocator()
            locator_pt.SetDataSet(self.Surface)
            locator_pt.BuildLocator()
            
            locator_ctr = vtk.vtkPointLocator()
            locator_ctr.SetDataSet(self.Centerlines)
            locator_ctr.BuildLocator()
            
            targetSeedIds = vtk.vtkIdList()
            outletSeedIds = vtk.vtkIdList()
            #print(self.BoundaryReference)
            for i in range(self.BoundaryReference.GetNumberOfPoints()):
                pt = self.BoundaryReference.GetPoint(i)
                br_vector = self.BoundaryReference.GetPointData().GetArray("BoundaryNormals").GetTuple(i)
                
                surf_ptId = locator_pt.FindClosestPoint(pt)
                ctr_ptId = locator_ctr.FindClosestPoint(pt)
                
                ctr_tangent = self.Centerlines.GetPointData().GetArray(self.FrenetTangentArrayName).GetTuple(ctr_ptId)
                
                direction = vtk.vtkMath.Dot(br_vector, ctr_tangent)
                if(direction < 0.0):
                    # outlet normal opposite centerline tangent
                    targetSeedIds.InsertNextId(surf_ptId)
                else:
                    outletSeedIds.InsertNextId(surf_ptId)
            
            mainBodySeedId = self.MainBodyPointId
        
        # set the point closest to the surface we want to keep
        mainBodyPoint = self.Surface.GetPoint(mainBodySeedId)            
        
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

        locator_ctr = vtk.vtkPointLocator()
        locator_ctr.SetDataSet(self.Centerlines)
        locator_ctr.BuildLocator()

        seeds_dict = { "inlets": targetSeedIds, "outlets" : outletSeedIds}

        ctr_points = vtk.vtkPoints()
        ctr_points_id = vtk.vtkUnsignedCharArray()
        ctr_points_id.SetNumberOfComponents(1)
        ctr_points_id.SetName("BoundaryRegion")

        ctr_pts_vector = vtk.vtkDoubleArray()
        ctr_pts_vector.SetName("BoundaryDirection")
        ctr_pts_vector.SetNumberOfComponents(3)

        vertices = vtk.vtkCellArray()
        #print(targetSeedIds.GetNumberOfIds(), outletSeedIds.GetNumberOfIds())
        seed_count = 0
        print(targetSeedIds.GetNumberOfIds(), outletSeedIds.GetNumberOfIds())

        for seed_type in seeds_dict.keys():
            #print(seeds_dict[seed_type].GetNumberOfIds())
            if (seeds_dict[seed_type].GetNumberOfIds() < 1):
                continue

            for i in range(seeds_dict[seed_type].GetNumberOfIds()):
                seedId = seeds_dict[seed_type].GetId(i)
                #print(seedId)

                #locator = vtk.vtkPointLocator()
                #locator.SetDataSet(clippedSurface)
                #locator.BuildLocator()

                seedPoint = self.Surface.GetPoint(seedId)
                #seedPointId = locator.FindClosestPoint(seedPoint)
                #surf_point = clippedSurface.GetPoint(seedPointId)

                centerlinePointId = copy.deepcopy(locator_ctr.FindClosestPoint(seedPoint))
                centerlinetangent = copy.deepcopy( self.Centerlines.GetPointData().GetArray(self.FrenetTangentArrayName).GetTuple(centerlinePointId))

                centerlinenormal = copy.deepcopy( self.Centerlines.GetPointData().GetArray(self.FrenetNormalArrayName).GetTuple(centerlinePointId))
                centerlinePoint = copy.deepcopy(self.Centerlines.GetPoint(centerlinePointId))

                ctr_points.InsertNextPoint(centerlinePoint)
                vertex = vtk.vtkVertex()
                vertex.GetPointIds().InsertId(0, seed_count)
                #vertex.GetPoints().SetPoint(seed_count, centerlinePoint)
                vertices.InsertNextCell(vertex)


                if( seed_type == "outlets"):
                    # outward facing normal
                    centerlinetangent = tuple( -p for p in centerlinetangent)
                    ctr_points_id.InsertNextTuple([int(1)])
                else:
                    ctr_points_id.InsertNextTuple([int(0)])
                ctr_pts_vector.InsertNextTuple(centerlinetangent) #inward pointing
                #planeEstimator = vtkvmtk.vtkvmtkPolyDataNormalPlaneEstimator()
                #planeEstimator.SetInputData(clippedSurface)
                #planeEstimator.SetOriginPointId(seedPointId)
                #planeEstimator.Update()

                plane = vtk.vtkPlane()
                plane.SetOrigin(centerlinePoint)
                plane.SetNormal(centerlinetangent)

                seamFilter = vtkvmtk.vtkvmtkTopologicalSeamFilter()
                seamFilter.SetInputData(clippedSurface)
                seamFilter.SetClosestPoint(seedPoint)
                seamFilter.SetSeamScalarsArrayName("SeamScalars")
                seamFilter.SetSeamFunction(plane)
                seamFilter.Update()

                #print("got here")
                # writer3 = vtk.vtkXMLPolyDataWriter()
                # writer3.SetFileName("seam_test_{0}.vtp".format(seed_count))
                # writer3.SetInputConnection(seamFilter.GetOutputPort())
                # writer3.Update()

                tree = vtk.vtkModifiedBSPTree()
                tree.SetDataSet(seamFilter.GetOutput())
                tree.BuildLocator()
                #intersect the locator with the line
                LineP0 = centerlinePoint
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
                    rotate.RotateWXYZ(theta,  centerlinetangent)
                    rotate.Update()

                    #print(dir(rotate))
                    #trans_m = vtk.vtkMatrix4x4()
                    #rotate.GetMatrix(trans_m)

                    out_vector = rotate.TransformVector(centerlinenormal)
                    LineP1 = [ c2 + count*dt*c1 for c2, c1 in zip(centerlinePoint, out_vector)]
                    #print(centerlinenormal, out_vector)
                    while ( code == 0 and count < 10000):
                        count += 1
                        code = tree.IntersectWithLine(LineP0, LineP1,
                                                      tolerance, IntersectPoints,
                                                      IntersectCells)
                        LineP1 = [ c2 + count*dt*c1 for c2, c1 in zip(centerlinePoint, out_vector)]
                    if(count > 10000 and code == 0):
                        print("no intersection")
                        continue

                    if (code != 0):
                        pt = IntersectPoints.GetPoint(0)
                        #pt = [ c2 + dt*c1 for c2, c1 in zip(pt, out_vector)] # add some buffer, may not need it
                        IntersectPointsList.InsertNextPoint(pt)
                        IntersectCellsList.InsertNextId(IntersectCells.GetId(0))
                        #print(IntersectPoints.GetPoint(0), IntersectCells.GetId(0) )

                # select = vtk.vtkImplicitSelectionLoop()
                # select.SetLoop(IntersectPointsList)
                # select.SetNormal(centerlinetangent)

                #set bool
                # Bool = vtk.vtkImplicitBoolean()
                # Bool.SetOperationTypeToIntersection() #SetOperationTypeToDifference()
                # Bool.AddFunction(select)
                # Bool.AddFunction(plane)

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

                # writer4 = vtk.vtkXMLPolyDataWriter()
                # writer4.SetFileName("selection_test_{0}.vtp".format(seed_count))
                # writer4.SetInputConnection(geometryFilter.GetOutputPort())
                # writer4.Update()

                # # grab all free edges and place them into a temporary polydata
                # num_sub_cells = geometryFilter.GetOutput().GetNumberOfCells()
                # neighbors = vtk.vtkIdList()
                # #neighbors.Allocate(vtk.VTK_CELL_SIZE)
                # free_edge_cells = set()
                # for cell_id in range(num_sub_cells):
                #     cellPointIds = vtk.vtkIdList()
                #     geometryFilter.GetOutput().GetCellPoints(cell_id, cellPointIds)
                #
                #     n_cell_pts = cellPointIds.GetNumberOfIds()
                #     for i in range(n_cell_pts):
                #         p1 = cellPointIds.GetId(i)
                #         p2 = cellPointIds.GetId((i+1) % n_cell_pts)
                #         geometryFilter.GetOutput().GetCellEdgeNeighbors(cell_id, p1, p2, neighbors)
                #         numNei = neighbors.GetNumberOfIds()
                #         if ( numNei < 1 ):
                #             free_edge_cells.add(cell_id)
                #             # newLines->InsertNextCell(2);
                #             # newLines->InsertCellPoint(p1);
                #             # newLines->InsertCellPoint(p2);
                #             # could add something to capture lines here, but currently
                #             # this is not yet necessary
                # print(len(free_edge_cells), num_sub_cells)

                ### try exluding to remesh subset, takes too long

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

                # writer4 = vtk.vtkXMLPolyDataWriter()
                # writer4.SetFileName("dist_test_{0}.vtp".format(seed_count))
                # writer4.SetInputConnection(geometryFilter.GetOutputPort())
                # writer4.Update()

                remeshsubset = vmtksubsetremesh()
                remeshsubset.CellEntityIdsArrayName = "RemeshCells"
                remeshsubset.ElementSizeMode = 'edgelengtharray'
                remeshsubset.TargetEdgeLengthArrayName = "PlaneDistance"
                remeshsubset.TargetEdgeLengthFactor = 0.8
                remeshsubset.NumberOfIterations = 10
                remeshsubset.Surface = geometryFilter.GetOutput()
                remeshsubset.ExcludeEntityIds = [0] #list(free_edge_cells)
                remeshsubset.Execute()

                # writer = vtk.vtkXMLPolyDataWriter()
                # writer.SetFileName("remesh_test_{0}.vtp".format(seed_count))
                # writer.SetInputData(remeshsubset.Surface)
                # writer.Update()

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

                # writer4 = vtk.vtkXMLPolyDataWriter()
                # writer4.SetFileName("remesh2_test_{0}.vtp".format(seed_count))
                # writer4.SetInputConnection(surfaceTriangulator2.GetOutputPort())
                # writer4.Update()

                # for cell_idx in hood_list:
                #     #print(cell_id)
                #     cellPointIds = vtk.vtkIdList()
                #     seamFilter.GetOutput().GetCellPoints(cell_idx, cellPointIds)
                #     for i in range(cellPointIds.GetNumberOfIds()):
                #         point_id = cellPointIds.GetId(i)
                #         pt = seamFilter.GetOutput().GetPoint(point_id)
                #         dist = plane.FunctionValue(pt)
                #
                #         if(dist > -1.0E-4 and dist < 1.0E-4 ):
                #             print("dist 1", dist)
                #
                #         if(dist > -1.0E-7 and dist < 1.0E-7 ):
                #             print("dist 2", dist)


                # locator = vtk.vtkPointLocator()
                # locator.SetDataSet(surfaceTriangulator2.GetOutput())
                # locator.BuildLocator()
                #
                # seedPoint =
                # seedPointId = locator.FindClosestPoint(seedPoint)
                # surf_point = surfaceTriangulator2.GetOutput().GetPoint(seedPointId)

                seamFilter2 = vtkvmtk.vtkvmtkTopologicalSeamFilter()
                seamFilter2.SetInputConnection(surfaceTriangulator2.GetOutputPort())
                seamFilter2.SetClosestPoint(IntersectPointsList.GetPoint(0))
                seamFilter2.SetSeamScalarsArrayName("SeamScalars")
                seamFilter2.SetSeamFunction(plane)
                seamFilter2.Update()

                # writer6 = vtk.vtkXMLPolyDataWriter()
                # writer6.SetFileName("seam_test2_{0}.vtp".format(seed_count))
                # writer6.SetInputConnection(seamFilter2.GetOutputPort())
                # writer6.Update()

                clipper = vtk.vtkClipPolyData()
                clipper.SetInputConnection(seamFilter2.GetOutputPort())
                #clipper.SetClipFunction(plane)
                clipper.GenerateClipScalarsOff()
                clipper.GenerateClippedOutputOn()
                #clipper.InsideOutOff()
                clipper.Update()

                writer7 = vtk.vtkXMLPolyDataWriter()
                writer7.SetFileName("clipper_test{0}.vtp".format(seed_count))
                writer7.SetInputConnection(clipper.GetOutputPort())
                writer7.Update()
                
                writer8 = vtk.vtkXMLPolyDataWriter()
                writer8.SetFileName("clipped_test{0}.vtp".format(seed_count))
                writer8.SetInputConnection(clipper.GetClippedOutputPort())
                writer8.Update()

                connectivity = vtk.vtkPolyDataConnectivityFilter()
                connectivity.SetInputConnection(clipper.GetOutputPort())
                connectivity.SetExtractionModeToClosestPointRegion()
                connectivity.SetClosestPoint(mainBodyPoint)
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
                seed_count += 1

                # writer7 = vtk.vtkXMLPolyDataWriter()
                # writer7.SetFileName("clipped_test{0}.vtp".format(seed_count))
                # writer7.SetInputData(clippedSurface)
                # writer7.Update()

        ctr_points_id.Squeeze()
        ctr_pts_vector.Squeeze()
        vertices.Squeeze()
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(ctr_points)
        polydata.SetVerts(vertices)
        polydata.GetCellData().AddArray(ctr_points_id)
        polydata.GetCellData().AddArray(ctr_pts_vector)
        #print(ctr_points.GetNumberOfPoints(), ctr_points_id.GetNumberOfTuples(), ctr_pts_vector.GetNumberOfTuples(), vertices.GetNumberOfCells())

        #writer7 = vtk.vtkXMLPolyDataWriter()
        #writer7.SetFileName("centerline_points.vtp")
        #writer7.SetDataModeToAscii()
        #writer7.SetInputData(polydata)
        #writer7.Update()
        
        self.BoundaryDirection = polydata

        self.Surface = clippedSurface

        if self.OwnRenderer:
            self.vmtkRenderer.Deallocate()


if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
