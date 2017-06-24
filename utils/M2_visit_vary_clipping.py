from visit import *
from visit_utils import *
import math
import csv, sys

# to run this script: visit -cli -nowin -s M2_visit_vary_clipping.py
name = 'M2_exp4_10km'
path = "/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx10km/Sigma_zstar/M2_exp4/VTK/"

RestoreSessionWithDifferentSources("M2_salt.session", 0, ("localhost:"+path+name+"-temp-*.vtk database","localhost:"+path+name+"-bathymetry.vtk","localhost:"+path+name+"-salt-*.vtk database"))

ResizeWindow(1, 800, 600)

# initial position
c0 = visit.View3DAttributes()
c0.viewNormal = (0.869777, 0.0963225, 0.483953)
c0.focus = (625, 1000, -300.038)
c0.viewUp =  (-0.475058, -0.101787, 0.874048)
c0.viewAngle = 30
c0.parallelScale = 1209.52
c0.nearPlane = -2430.67
c0.farPlane = 2430.67
c0.imagePan = (-0.0194483, -0.0190144)
c0.imageZoom = 1.35171
c0.perspective = 0
c0.eyeAngle = 2
c0.centerOfRotationSet = 0
c0.centerOfRotation = (625, 1000, -300.038)
c0.axis3DScaleFlag = 0
c0.axis3DScales = (1, 1, 1)
c0.shear = (0, 0, 1)

# set view to these values
SetView3D(c0)

# number of snapshots
n=TimeSliderGetNStates()
# Initial Clipping AAttributes
pos0 = 1250 # initial pos
ClipAtts = visit.ClipAttributes()
ClipAtts.plane1Status = 1
ClipAtts.plane1Origin = (pos0, 0, 0)
ClipAtts.plane1Normal = (1, 0, 0)
ClipAtts.planeInverse = 0
ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1  # None, Plane1, Plane2, Plane3
ClipAtts.center = (0, 0, 0)
ClipAtts.radius = 1
ClipAtts.sphereInverse = 0
SetOperatorOptions(ClipAtts, 2)

# Set basic save options
s = SaveWindowAttributes()
#
# The 'family' option controls if visit automatically adds a frame number to
# the rendered files. For this example we will explicitly manage the output name.
#
s.family = 0
s.format = s.PNG
s.screenCapture = 0
s.quality = 100
s.progressive = 0
s.width = 900
s.height = 600


# indices controlling clip
ti = 200; tf = 230; trho = 300
pos1 = 600 # final clip
cstep = (pos0-pos1)/(tf-ti) # clip step
pos = pos0 # initial pos
for i in range(n):
    if i>=ti and i<tf:
       pos = pos - cstep
       print 'pos0-cstep = ', pos
       #Clip
       ClipAtts.plane1Origin = ((pos), 0, 0)
    elif i>=tf:
       ClipAtts.plane1Origin = (pos1, 0, 0)

    if i == trho:
       #SetActivePlots((2, 5))
       #SetActivePlots(5)
       SetActivePlots((1, 3))
       SetActivePlots(3)
       HideActivePlots()
    elif i<trho:
       ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1
       ClipAtts.plane1Status = 1
       ClipAtts.plane2Status = 0
       # update visit
       SetOperatorOptions(ClipAtts, 0)

    #before we render the result, explicitly set the filename for this render
    s.fileName = name+"-%04d.png" % i
    SetSaveWindowAttributes(s)
    # render the image to a PNG file
    SetTimeSliderState(i) # advance timeslider
    DrawPlots()
    SaveWindow()

# show time for last frame
#banner.text = str("Time = %5.2f days" % (time0-dt))
View3DAttributes()
sys.exit()
