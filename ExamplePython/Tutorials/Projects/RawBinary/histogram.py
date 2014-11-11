try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

ProcessIdScalars1 = GetActiveSource()
Histogram1 = Histogram()

Histogram1.CustomBinRanges = [-0.97492790222167969, 0.97492790222167969]
Histogram1.SelectInputArray = ['POINTS', 'Normals']

Histogram1.SelectInputArray = ['POINTS', 'random number']

AnimationScene1 = GetAnimationScene()
RenderView1 = GetRenderView()
XYBarChartView1 = CreateBarChartView()
XYBarChartView1.ViewTime = 0.0

DataRepresentation1 = Show()
DataRepresentation1.XArrayName = 'bin_extents'
DataRepresentation1.SeriesVisibility = ['bin_extents', '0', 'Normals_total (0)', '0', 'Normals_total (1)', '0', 'Normals_total (2)', '0', 'Normals_total (Magnitude)', '0', 'Normals_average (0)', '0', 'Normals_average (1)', '0', 'Normals_average (2)', '0', 'Normals_average (Magnitude)', '0', 'ProcessId_total', '0', 'ProcessId_average', '0', 'vtkOriginalIndices', '0', 'bin_values', '1']
DataRepresentation1.AttributeType = 'Row Data'
DataRepresentation1.UseIndexForXAxis = 0

AnimationScene1.ViewModules = [ RenderView1, XYBarChartView1 ]

Render()
