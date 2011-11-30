try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

c02_200_127_gdg = GadgetReader( AaaFileName='/Users/corbett/Local/BigData/HalosGoingNotts/aquarius/A-5/c02_200_127.gdg' )

c02_200_127_gdg.PointArrays = ['Energy', 'Mass', 'Type', 'Velocity']

c02_200_127_gdg.LUnit = 9.9999999999999995e-07
c02_200_127_gdg.FilePrefix = '/Users/corbett/Local/BigData/HalosGoingNotts/aquarius/A-5/snap_C02_200_127'
c02_200_127_gdg.Swap = 0
c02_200_127_gdg.PointArrays = ['Energy', 'Mass', 'Type', 'Velocity']

RenderView1 = GetRenderView()
DataRepresentation3 = Show()
DataRepresentation3.EdgeColor = [0.0, 0.0, 0.50000762951094835]

RenderView1.CameraPosition = [-1.1375389426208039e+35, -3.5254681548971297e+36, 2.221829504886341e+39]
RenderView1.CameraClippingRange = [1.5586970978390507e+39, 3.0867182996554522e+39]
RenderView1.CameraFocalPoint = [-1.1375389426208039e+35, -3.5254681548971297e+36, -1.0959873841862829e+37]
RenderView1.CameraParallelScale = 5.7788841491748422e+38
RenderView1.CenterOfRotation = [-1.1375389426208039e+35, -3.5254681548971297e+36, -1.0959873841862829e+37]

Render()
