try: paraview.simple
except: from paraview.simple import *
ROOT_DIR='/Users/corbett/Documents/Projects/EuroPython/'
def main():
  progSource =\
    ProgrammableSource(Script=\
    open(ROOT_DIR+'/Tutorials/Projects/CustomReader/SimplePoints.py').read())
  progSource.UpdatePipeline()
  glyph = Glyph(GlyphType='Sphere',Scalars = ['POINTS', 'Temperature'],     ScaleMode = 'off')
  dr=Show()
  dr.ColorArrayName='Temperature'
  dr.LookupTable=GetLookupTableForArray('Temperature', 1 )
  Render()
  ResetCamera()
  WriteImage(ROOT_DIR+'/Tutorials/Screenshots/CustomReader.png')

if __name__ == '__main__':
  main()
  
  
  
  

