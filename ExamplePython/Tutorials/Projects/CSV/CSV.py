try: paraview.simple
except: from paraview.simple import *
import random
ROOT_DIR='/Users/corbett/Documents/Projects/EuroPython/'

def gen_simple_csv():
  f=open(ROOT_DIR+'Tutorials/Data/simple.csv','w')
  f.write('x,y,z,random\n')
  for x in range(10):
    for y in range(10):
      for z in range(10):
        f.write('%d,%d,%d,%f\n' % (x,y,z,random.uniform(0,10)))


def main():
  # Reading in from file 
  points = CSVReader(FileName=[ROOT_DIR+'/Tutorials/Data/simple.csv'])
  points.UpdatePipeline()
  # Making this spreadsheet into coordinates we can visualize
  table = TableToPoints()
  table.XColumn = 'x'
  table.YColumn = 'y'
  table.ZColumn = 'z'
  table.UpdatePipeline()
  # Instead of points lets visualize things as a sphere
  sphereglyph = Glyph(GlyphType='Sphere')
  dr=Show()
  dr.ColorArrayName='random'
  Render()
  ResetCamera()
  #Now Let's Save the Image into a Png
  WriteImage(ROOT_DIR+'/Tutorials/Screenshots/CSV.png')
  
  #Now let's see what process these are on
  pids = ProcessIdScalars()
  dr=Show()
  dr.ColorArrayName='ProcessId'
  Render()
  #Now Let's Save the Image into a Png
  WriteImage(ROOT_DIR+'/Tutorials/Screenshots/CSV-processorIDs-before-D3.png')

  # Uh-oh, everything is of the same color, i.e. on the same process! This is because the comma-separated value filter takes in an ascii file and is a fundamentally serial reader, so all the data is currently on process 0. We can use one of ParaView's filters called the Data Distribution Filter or D3 for short to change this
  d3= D3()
  dr=Show()

  #Finally let's visualize the process ids this time around
  pids = ProcessIdScalars()
  dr=Show()
  dr.ColorArrayName='ProcessId'
  Render()
  WriteImage(ROOT_DIR+'/Tutorials/Screenshots/CSV-processorIDs-after-D3.png')  
if __name__ == '__main__':
  main()