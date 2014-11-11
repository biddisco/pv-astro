from paraview import vtk 
import os 
# Abstractions
class Array:
  def __init__(self,name,values):
    self.vtkarray=vtk.vtkDoubleArray()
    self.vtkarray.SetName(name)
    try:
      self.vtkarray.SetNumberOfComponents(len(values[0]))
    except:
      self.vtkarray.SetNumberOfComponents(1)
    for val in values:
      self.addValue(val)
  #def __init__(self,name): # can't do this in python
  #  self.__init__(name,[])
  def addValue(self,val):
    self.array.InsertNextValue(val)
  def addTo(self,output):
      output.GetPointData().AddArray(self.vtkarray)
class DataSet:
  def __init__(self,points,arrays):
    self.arrays=arrays
    self.points=vtk.vtkPoints()
    for pnt in points:
      self.addPoint(pnt[0],pnt[1],pnt[2])
  #def __init__(self,arrays): #can't do this in python
  #  self.__init__([],arrays)
  def addPoint(self,x,y,z):
    self.points.InsertNextPoint(x,y,z)
  def addValueToArray(self,array_name,value):
    self.arrays[array_name].addValue(value)
  def addTo(self,output):
    output.SetPoints(self.points)
    for ary in self.arrays:
      ary.addTo(output)
####
# Actual work
####
# Our hardcoded dataset
points_to_create=[(0,0,0),(1,1,1),(2,2,2),(3,3,3)]
array=Array('Temperature',[0.0,1.1,2.2,3.3])
data=DataSet(points_to_create,[array])

# This is where we'll place our dataset
output = self.GetOutput()
data.addTo(output)

