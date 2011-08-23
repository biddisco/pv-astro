from math import pi
ramses_file=GetActiveSource()

# CGS units needed
kpc_in_cm=3.08568025*10**21 #cgs_units
pc_in_cm=3.08568025*10**18
km_in_cm=10**5
Gyr_in_s=3.1536*10**16
yr_in_s=3.1553*10**7
msol_in_g=1.98892*10**33
# other important constants
G=4.3*10**-3 #pc/msol*(km/s)^2
h0=73 #km/s/Mpc, default could change

boxlen=1 #default
unit_l=1 #default-should change!                                                                                                                  
unit_m=1 #default-should change!
unit_d=1 #default-should change!
unit_t=1 #default-should change!                                                                                                            
for i in range(ramses_file.FieldData.NumberOfArrays):
    arr=ramses_file.FieldData.GetArray(i)
    if arr.Name=='unit_l':
        unit_l=arr.GetRange()[0]
        print 'unit_l=%f'%unit_l
    elif arr.Name=='unit_m':
        unit_m=arr.GetRange()[0]
        print 'unit_m=%f'%unit_m
    elif arr.Name=='unit_d':
        unit_d=arr.GetRange()[0]
        print 'unit_d=%e'%unit_d
    elif arr.Name=='unit_t':
        unit_t=arr.GetRange()[0]
        print 'unit_t=%f'%unit_t
    elif arr.Name=='boxlen':
        boxlen=arr.GetRange()[0]
        print 'boxlen=%f'%boxlen
    elif arr.Name=='h0':
        h0=arr.GetRange()[0]
        print 'h0=%f'%boxlen
conv_factor=boxlen*unit_l/kpc_in_cm
print 'conv_factor=%f' % conv_factor
rho_crit=3.*(h0/10**3)**2/(8*pi*G/10**3) #msol/kpc^3    
print '200*rho_crit=%f' % (200*rho_crit)                                                                            

print 'Converting coordinates from simulation units to kpc'
calc = Calculator(AttributeMode = 'point_data',CoordinateResults = 1,Function ='%f*coords' % conv_factor)
calc.UpdatePipeline()
Show()
Render()

print 'Converting velocities from simulation units to km/s'
calc = Calculator(AttributeMode = 'point_data',
  ResultArrayName = 'velocitykmpers',
  Function ='velocity*(%f)/(%f)/(%f)' % (unit_l,unit_t,km_in_cm))
calc.UpdatePipeline()
Show()
Render()

print 'Converting masses from simulation units to Msol'
calc = Calculator(AttributeMode = 'point_data',
  ResultArrayName = 'massmsol',
  Function ='mass*(%f)^3*(%e)/(%f)' % (unit_l,unit_d,msol_in_g))
calc.UpdatePipeline()
Show()
Render()

print 'Converting ages to Gyr (useful for average stellar ages)'
calc = Calculator(AttributeMode = 'point_data',
  ResultArrayName = 'ageinGyr',
  Function ='-1*age*(%f)/(%f)' % (unit_t,Gyr_in_s))
calc.UpdatePipeline()
Show()
Render()

