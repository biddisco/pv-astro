from paraview import vtk
# where we'll put our data
output=self.GetOutput()

# Our data: points and arrays
# Points
points=vtk.vtkPoints()
for point in [(0,0,0),(1,1,1),(2,2,2),(3,3,3)]:
  points.InsertNextPoint(point[0],point[1],point[2])

# One Scalar Array
temp=vtk.vtkDoubleArray()
temp.SetName('Temperature')
temp.SetNumberOfComponents(1)

for val in [0.0,1.1,2.2,3.3]:
  temp.InsertNextValue(val)

# Placing these in the output dataset
output.SetPoints(points)
output.GetPointData().AddArray(temp)
# There we are!
  

  
