from math import pi
ROOT_DIR='/Users/corbett/Documents/Projects/pvaddons/'
# Reading in file
ramses_file = RamsesReader(AaaFileName=\
    '/home/itp/teyssier/lustre/AqC6/output_00193/info_00193.txt' ,
   PointArrays =\
    ['Age', 'Eps', 'Hsmooth', 'Mass', 'Metals', 'Potential', 
    'Rho', 'Temperature', 'Tform', 'Type', 'Velocity'],
   ParticleMassGuess = 3000000000.0) # particle mass guess large as we don't want to bother creating dark particles
ramses_file.Show()
Render()

# Converting units to units we are happy with
unitConversionCode=open(ROOT_DIR+'ParaViz/PythonScripts/unit_conversion.py').read()
exec unitConversionCode

# Selecting only dark particles
thresh=Threshold(Scalars = ['POINTS', 'type'],ThresholdRange = [0.0, 0.0])
Show()
Render()

# Now converting
center=[0.567804, 0.586049, 0.559157] #hardcoded, got from PV-meshless which can in principle be scripted
#center=[77643.001980356319, 80137.870933607104, 76460.588615666857]#=map(lambda x: conv_factor*x,center)
if not center:
  ## TODO: here use pv-meshless to script and calculate the center=
  center=[0.0,0.0,0.0]
  ppf = ParticlePartitionFilter()
  sph = SPHProbeCustomsource( Probe=ppf,
    ComputeDensityFromNeighbourVolume = 1,
    MassScalarArray = 'mass',
    HsmoothinglengthArray = 'eps',
    DensityScalarArray = 'Not available',
    VolumeScalarArray = 'Not available')
  # Extract the density out of the result
  
  # Threshold by this density
  
  # Select the coordinates of the particle with this density
  
  # This is our center, set the center to this






# Selecting the inner 20 kpc for further analysis
clip = Clip( ClipType="Sphere",InsideOut = 1 )
clip.ClipType.Center = map(lambda x: conv_factor*x,center)
clip.ClipType.Radius = 20.0
Show()
Render()
ResetCamera()

# Creating a Glyph
glyph = Glyph(RandomMode = 0,
  GlyphType = "2D Glyph"
  ScaleMode = 'off'
  MaskPoints = 0)
glyph.GlyphType.GlyphType = 'Vertex'
Show()
Render()
# Doing a profile

profile = Profile( ProbeType="Fixed Radius Point Source",
    SelectInputArray = ['POINTS', 'mass msol'])
profile.ProbeType.Center = map(lambda x: conv_factor*x,center)

# Computing the virial radius
#p_crit=