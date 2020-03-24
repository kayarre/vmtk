#!/usr/bin/env python

## Modification by Peter Patalano  to the below software
## Modification Summary
## * Added AutoClip functionality which allows clips from a mesh to be stored then used again on a future mesh
## * Useful for meshes that have undergone a stenosis following CFD and need to remove old extensions to place new ones
## *
## *
## Packages Modified
## * vmtksurfaceclipper - see above
## * vmtksurfacereader - added MultiBlockDataSet functionality
## * vmtksurfacewriter - added MultiBlockDataSet functionality
## *




## Program:   VMTK
## Module:    $RCSfile: vmtksurfaceclipper.py,v $
## Language:  Python
## Date:      $Date: 2006/05/26 12:35:13 $
## Version:   $Revision: 1.9 $

##   Copyright (c) Luca Antiga, David Steinman. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

from __future__ import absolute_import #NEEDS TO STAY AS TOP LEVEL MODULE FOR Py2-3 COMPATIBILITY
import vtk
import sys

# Debugging
# import pdb

from vmtk import vmtkrenderer
from vmtk import pypes

class vmtkSurfaceClipper(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.Surface = None
        self.ClippedSurface = None
        self.CutLines = None
        self.vmtkRenderer = None
        self.OwnRenderer = 0

        self.WidgetType = 'box'

        self.Actor = None
        self.ClipWidget = None
        self.ClipFunction = None
        self.CleanOutput = 1
        self.Transform = None

        self.InsideOut = 0

        self.Interactive = 1
        self.ClipArrayName = None
        self.ClipValue = 0.0
        
        # Globals for our added functionality
        self.ClipsIn = None 
        self.ClipsOut = None
        self.ClipCollection = None

        self.SetScriptName('petersurfaceclipper')
        self.SetScriptDoc('interactively clip a surface with a box')
        self.SetInputMembers([
            ['Surface','i','vtkPolyData',1,'','the input surface','vmtksurfacereader'],
            # The new input member which will be a MultiBlock file
            ['ClipsIn','iclip','vtkMultiBlockDataSet',1,'','the input clips file','vmtksurfacereader'],
            ['WidgetType','type','str',1,'["box","sphere"]','type of widget used for clipping'],
            ['Transform','transform','vtkTransform',1,'','the widget transform, useful in case of piping of multiple clipping scripts'],
            ['CleanOutput','cleanoutput','bool',1,'','toggle cleaning the unused points'],
            ['InsideOut','insideout','bool',1,'','toggle switching output and clipped surfaces'],
            ['Interactive','interactive','bool',1,'','toggle clipping with a widget or with a pre-defined scalar field'],
            ['ClipArrayName','array','str',1,'','name of the array with which to clip in case of non-interactive clipping'],
            ['ClipValue','value','float',1,'','scalar value at which to perform clipping in case of non-interactive clipping'],
            ['vmtkRenderer','renderer','vmtkRenderer',1,'','external renderer']
            ])
        self.SetOutputMembers([
            ['Surface','o','vtkPolyData',1,'','the output surface','vmtksurfacewriter'],
            ['ClippedSurface','oclipped','vtkPolyData',1,'','the clipped surface','vmtksurfacewriter'],
            ['CutLines','ocutlines','vtkPolyData',1,'','the cutlines','vmtksurfacewriter'],
            # The new otput member which will be a MultiBlock file
            ['ClipsOut','oclip','vtkMultiBlockDataSet',-1,'','the output surfaces associated with the clips made','vmtksurfacewriter'],
            ['Transform','otransform','vtkTransform',1,'','the output widget transform']
            ])
    
    def AutoClip(self):
        # Check to see if ClipsIn is provided
        if self.ClipsIn == None:
            self.PrintError('Error: no AutoClip without ClipsIn')
        
        
        import pdb; pdb.set_trace() 
        # Init AutoClipper/Function
        AutoClipper = vtk.vtkClipPolyData()
        AutoClipFunction = vtk.vtkPlanes()

        AutoClipper.SetInputData(self.Surface)
        AutoClipper.GenerateClippedOutputOn()
        AutoClipper.SetInsideOut(self.InsideOut)
        AutoClipper.SetClipFunction(AutoClipFunction)
        
        # Perform ClipFunction while ClipsIn is not empty
        while not self.ClipsIn.GetNumberOfBlocks() == 0:
            Bounds = self.NextBlockBounds()
            self.PrintLog('Bounds for AutoClip are: ' + str(Bounds))
            AutoClipFunction.SetBounds(Bounds)
            AutoClipper.Update()
            self.Surface.DeepCopy(AutoClipper.GetOutput())
        
        self.PrintLog('AutoClip complete.')

    def NextBlockBounds(self):
        if self.ClipsIn == None:
            self.PrintError('Error: no clip input provided')
        
        # Get the index for the last block in the dataset
        LastBlock = self.ClipsIn.GetNumberOfBlocks() - 1
        
        # Get data for the last block
        BlockData = vtk.vtkPolyData()
        BlockData.DeepCopy(self.ClipsIn.GetBlock(LastBlock))
        self.ClipsIn.RemoveBlock(LastBlock)
        Bounds = BlockData.GetBounds()
        return Bounds
        
    def GetClipData(self):
        #Initialize ClipCollection
        if self.ClipCollection == None:
            self.ClipCollection = vtk.vtkPolyDataCollection()

        #Get PolyData of the current box, requires deep copy
        dummy = vtk.vtkPolyData()
        ClipPolyData = vtk.vtkPolyData()
        self.ClipWidget.GetPolyData(dummy)
        ClipPolyData.DeepCopy(dummy)
        
        #Add ClipPolyData to the ClipCollection
        self.ClipCollection.AddItem(ClipPolyData)
        
        #Display a running list of clips made
        self.ClipCollection.InitTraversal()
        self.PrintLog('All bounds are: ')
        for i in range(self.ClipCollection.GetNumberOfItems()):
            a = self.ClipCollection.GetNextItem()
            self.PrintLog(str(a.GetBounds()))

    def ClipCallback(self, obj):
        if self.ClipWidget.GetEnabled() != 1:
            return
        
        if self.WidgetType == "box":
            #self.ClipWidget.SetPlaceFactor(1)
            #self.ClipWidget.PlaceWidget(self.NextBlockBounds())
            self.Display()
            self.ClipWidget.GetPlanes(self.ClipFunction)

        elif self.WidgetType == "sphere":
            self.ClipWidget.GetSphere(self.ClipFunction)
        self.GetClipData()
        self.Clipper.Update()
        self.Surface.DeepCopy(self.Clipper.GetOutput())
        self.ClippedSurface.DeepCopy(self.Clipper.GetClippedOutput())
        self.Cutter.Update()
        self.CutLines.DeepCopy(self.Cutter.GetOutput())
        self.ClipWidget.Off()

    def InteractCallback(self, obj):
        pass
    #    if self.ClipWidget.GetEnabled() == 1:
    #        self.ClipWidget.SetEnabled(0)
    #    else:
    #        self.ClipWidget.SetEnabled(1)

    def Display(self):
        self.ClipWidget.SetInputData(self.Surface)
        self.ClipWidget.PlaceWidget()
        
        # If ClipsIn provided, display given bounds of last block
        #if self.ClipsIn and self.WidgetType == "box":
        #    self.ClipWidget.SetPlaceFactor(1)
        #    self.ClipWidget.PlaceWidget(self.NextBlockBounds())
        #    self.ClipWidget.On()
        
        if self.Transform and self.WidgetType == "box":
            self.PrintLog('ClipsIn is None')
            self.ClipWidget.SetTransform(self.Transform)
            self.ClipWidget.On()

        #self.vmtkRenderer.RenderWindowInteractor.Initialize()
        self.vmtkRenderer.Render()
        self.PrintLog('End of Display()')
       
        #self.vmtkRenderer.RenderWindowInteractor.Start()

    def Execute(self):

        if self.Surface == None:
            self.PrintError('Error: no Surface to for Execute().')

        self.Clipper = vtk.vtkClipPolyData()
        self.Clipper.SetInputData(self.Surface)
        self.Clipper.GenerateClippedOutputOn()
        self.Clipper.SetInsideOut(self.InsideOut)

        if self.Interactive:

            if self.WidgetType == "box":
                self.ClipFunction = vtk.vtkPlanes()
            elif self.WidgetType == "sphere":
                self.ClipFunction = vtk.vtkSphere()
       
            self.Clipper.SetClipFunction(self.ClipFunction)

            self.Cutter = vtk.vtkCutter()
            self.Cutter.SetInputData(self.Surface)
            self.Cutter.SetCutFunction(self.ClipFunction)

            self.ClippedSurface = vtk.vtkPolyData()
            self.CutLines = vtk.vtkPolyData()

            if not self.vmtkRenderer:
                self.vmtkRenderer = vmtkrenderer.vmtkRenderer()
                self.vmtkRenderer.Initialize()
                self.OwnRenderer = 1

            self.vmtkRenderer.RegisterScript(self) 

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(self.Surface)
            mapper.ScalarVisibilityOff()
            self.Actor = vtk.vtkActor()
            self.Actor.SetMapper(mapper)
            self.vmtkRenderer.Renderer.AddActor(self.Actor)
            if self.WidgetType == "box":
                self.ClipWidget = vtk.vtkBoxWidget()
                self.ClipWidget.GetFaceProperty().SetColor(0.6,0.6,0.2)
                self.ClipWidget.GetFaceProperty().SetOpacity(0.25)
            elif self.WidgetType == "sphere":
                self.ClipWidget = vtk.vtkSphereWidget()
                self.ClipWidget.GetSphereProperty().SetColor(0.6,0.6,0.2)
                self.ClipWidget.GetSphereProperty().SetOpacity(0.25)
                self.ClipWidget.GetSelectedSphereProperty().SetColor(0.6,0.0,0.0)
                self.ClipWidget.GetSelectedSphereProperty().SetOpacity(0.75)
                self.ClipWidget.SetRepresentationToSurface()
                self.ClipWidget.SetPhiResolution(20)
                self.ClipWidget.SetThetaResolution(20)

            self.ClipWidget.SetInteractor(self.vmtkRenderer.RenderWindowInteractor)
            # Try AutoClip here
            if not self.ClipsIn == None:
                self.AutoClip()

            self.vmtkRenderer.AddKeyBinding('space','Clip.',self.ClipCallback)
            self.vmtkRenderer.AddKeyBinding('i','Interact.',self.InteractCallback)
            self.Display()

            self.Transform = vtk.vtkTransform()
            if self.WidgetType == "box":
                self.ClipWidget.GetTransform(self.Transform)

            # Add items from ClipCollection to MultiBlockDataSet
            self.ClipsOut = vtk.vtkMultiBlockDataSet()
            ClipTotal = self.ClipCollection.GetNumberOfItems() 
            BlockNum = 0
            self.ClipCollection.InitTraversal()
            while BlockNum < ClipTotal:
                blockin = vtk.vtkPolyData()
                blockin.DeepCopy(self.ClipCollection.GetNextItem())
                self.PrintLog('The bounds of blockin are: ' + str(blockin.GetBounds()))
                self.ClipsOut.SetBlock(BlockNum, blockin)
                BlockNum += 1
                

            if self.OwnRenderer:
                self.vmtkRenderer.Deallocate()

        else:

            self.Surface.GetPointData().SetActiveScalars(self.ClipArrayName)

            self.Clipper.GenerateClipScalarsOff()
            self.Clipper.SetValue(self.ClipValue)
            self.Clipper.Update()

            self.Cutter = vtk.vtkContourFilter()
            self.Cutter.SetInputData(self.Surface)
            self.Cutter.SetValue(0,self.ClipValue)
            self.Cutter.Update()

            self.Surface = self.Clipper.GetOutput()
            self.ClippedSurface = self.Clipper.GetClippedOutput()

            self.CutLines = self.Cutter.GetOutput()

        if self.CleanOutput == 1:

            cleaner = vtk.vtkCleanPolyData()
            cleaner.SetInputData(self.Surface)
            cleaner.Update()
            self.Surface = cleaner.GetOutput()

            cleaner = vtk.vtkCleanPolyData()
            cleaner.SetInputData(self.ClippedSurface)
            cleaner.Update()
            self.ClippedSurface = cleaner.GetOutput()

            cleaner = vtk.vtkCleanPolyData()
            cleaner.SetInputData(self.CutLines)
            cleaner.Update()
            stripper = vtk.vtkStripper()
            stripper.SetInputConnection(cleaner.GetOutputPort())
            stripper.Update()
            self.CutLines = stripper.GetOutput()



if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
