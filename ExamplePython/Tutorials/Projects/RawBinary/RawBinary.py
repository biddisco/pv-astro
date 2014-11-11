import struct
import random
try: paraview.simple
except: from paraview.simple import *
ROOT_DIR='/Users/corbett/Documents/Projects/EuroPython/'
class BinaryWriter:
  def __init__(self,binfile,datatype):
      """
      binfile: path to file to write to
      datatype: one of format strings in
        http://docs.python.org/release/2.5.2/lib/module-struct.html
      """
      self.binfile=open(binfile, 'wb')
      self.datatype=datatype 
  def write(self,data):
    for datum in data:
      datum=struct.pack(self.datatype,datum)
      self.binfile.write(datum)

def gen_simple_raw():
  writer=BinaryWriter('coordinates.raw','i')
  for z in range(10):    
    for y in range(10):
      for x in range(10):
        writer.write([x,y,z])

def gen_random_number_raw():
  writer=BinaryWriter('random.raw','d')
  for z in range(10):    
    for y in range(10):
      for x in range(10):
        writer.write([random.random()])

def vis_bin_file():
  random_raw = ImageReader(FilePrefix=ROOT_DIR+'/Tutorials/Data/random.raw',
    ScalarArrayName = 'random number',
    DataExtent = [0, 9, 0, 9, 0, 9],
    DataByteOrder = 'LittleEndian',
    DataScalarType = 'double')
  random_raw.UpdatePipeline()

  #Structured data such as this doesn't automatically get a display in the GUI, we have to create a glyph.
  glyph = Glyph( GlyphType="Sphere", 
     Scalars = ['POINTS', 'random number'],
     ScaleMode = 'scalar'
     )
  dr=Show()
  dr.ColorArrayName='random number'
  Render()
  WriteImage(ROOT_DIR+'/Tutorials/Screenshots/RawBinary.png')
  # Now if we want to select some subset of data
  # say the random numbers between 0.75 and 1
  thresh = Threshold(Scalars = ['POINTS', 'random number'],
    ThresholdRange = [.75, 1.0])
  dr=Show()
  dr.ColorArrayName='random number'
  Render()
  ResetCamera()
  # Saving a screenshot
  WriteImage(ROOT_DIR+'/Tutorials/Screenshots/RawBinary-topquart.png')
  
  # Getting information about the number of points after our threshold
  di=thresh.GetDataInformation()
  print 'data set has %d points' % (di.DataInformation.GetNumberOfPoints())
  
  #Finally let's visualize the process ids this time around
  pids = ProcessIdScalars()
  dr=Show()
  dr.ColorArrayName='ProcessId'
  Render()
  
  
  
  
     
if __name__ == '__main__':
  #gen_simple_raw()
  #gen_random_number_raw()
  vis_bin_file()