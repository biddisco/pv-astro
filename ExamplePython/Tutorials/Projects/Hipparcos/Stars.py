try: paraview.simple
except: from paraview.simple import *
from math import pi
ROOT_DIR='/Users/corbett/Documents/Projects/EuroPython/'
def main():
  hipparcos = \
    CSVReader( FileName=[ROOT_DIR+'/Tutorials/Data/HipparcosPrune.csv'])
  # Dummy table to points filter so that we can
  # use filters that work only on point data, such as calc
  table=TableToPoints(
    XColumn = 'StarID',
    YColumn = 'StarID',
    ZColumn = 'StarID'
  )
  table.UpdatePipeline()
  # This will give us x,y,z from ascencion and declination
  calc=Calculator(
    Function ='cos(Dec)*cos(RA)/%f*iHat+cos(Dec)*sin(RA)/%f*jHat+sin(Dec)/%f*kHat'%(pi,pi,pi),
    CoordinateResults = 1
  )
  calc.UpdatePipeline()
  glyph = Glyph(GlyphType='2D Glyph',Scalars = ['POINTS', 'AbsMag'])
  glyph.GlyphType.GlyphType = 'Vertex'
  dr=Show()
  dr.ColorArrayName = 'AbsMag'
  dr.LookupTable=GetLookupTableForArray('AbsMag', 1 )

  Render()
  ResetCamera()
  WriteImage(ROOT_DIR+'/Tutorials/Screenshots/Hipparcos.png')

if __name__ == '__main__':
  main()