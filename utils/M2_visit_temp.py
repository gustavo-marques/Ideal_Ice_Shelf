from visit import *
from visit_utils import *
import math, sys
import csv

#run this script as follows: visit -cli -nowin -s M2_visit_temp.py

RestoreSessionWithDifferentSources("M2_temp_cliped.session", 0, ("localhost:/ptmp/gmm/M2_exp0_2km-*.vtk database","localhost:/archive/gmm/Ideal_ice_shelf/Mode2/dx2km/Sigma_zstar/M2_exp0/VTK/M2_exp0_2km-bathymetry.vtk","localhost:/archive/gmm/Ideal_ice_shelf/Mode2/dx2km/Sigma_zstar/M2_exp0/VTK/M2_exp0_2km-ice-shelf.vtk"))
ResizeWindow(1, 800, 600)
# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0.89164, 0.0693782, 0.447398)
View3DAtts.focus = (625, 1000, -299.888)
View3DAtts.viewUp = (-0.436022, -0.134512, 0.889827)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 1209.52
View3DAtts.nearPlane = -2427.75
View3DAtts.farPlane = 2427.75
View3DAtts.imagePan = (-0.0194483, -0.0190144)
View3DAtts.imageZoom = 1.35171
View3DAtts.perspective = 0
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (625, 1000, -300.038)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state

# Set basic save options
s = SaveWindowAttributes()
#
# The 'family' option controls if visit automatically adds a frame number to 
# the rendered files. For this example we will explicitly manage the output name.
#
s.family = 0
s.format = s.PNG 
# set the width of the output image
s.screenCapture = 0
s.quality = 100
s.progressive = 0
s.width = 900
s.height = 600

# number of snapshots
n=TimeSliderGetNStates()

# Initial Clipping AAttributes
pos0 = 1250 # initial pos
#ClipAtts = visit.ClipAttributes()
ClipAtts = GetOperatorOptions(2)
#ClipAtts.quality = visit.ClipAtts.Fast  # Fast, Accurate
#ClipAtts.funcType = visit.ClipAtts.Plane  # Plane, Sphere
ClipAtts.plane1Status = 1
ClipAtts.plane1Origin = (pos0, 0, 0)
ClipAtts.plane1Normal = (1, 0, 0)
ClipAtts.planeInverse = 0
ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1  # None, Plane1, Plane2, Plane3
ClipAtts.center = (0, 0, 0)
ClipAtts.radius = 1
ClipAtts.sphereInverse = 0
#SetOperatorOptions(ClipAtts, 2)

# indices controlling clip
ti = 200; tf = 230; trho = 300
pos1 = 500 # final clip
cstep = (pos0-pos1)/(tf-ti) # clip step
pos = pos0 # initial pos
for i in range(n):
    #if i>=ti and i<tf:
    #   pos = pos - cstep
    #   print 'pos0-cstep = ', pos
    #   #Clip
    #   ClipAtts.plane1Origin = ((pos), 0, 0)
    #elif i>=tf:
    #   ClipAtts.plane1Origin = (pos1, 0, 0)
    #
    #if i == trho:
    #   SetActivePlots((2, 5))
    #   SetActivePlots(5)
    #   HideActivePlots()
    #elif i<trho:
    #   ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1
    #   ClipAtts.plane1Status = 1
    #   ClipAtts.plane2Status = 0
    #   # update visit
    #   visit.SetOperatorOptions(ClipAtts, 2)

    #before we render the result, explicitly set the filename for this render
    s.fileName = "example_%04d.png" % i
    SetSaveWindowAttributes(s)
    # render the image to a PNG file
    SetTimeSliderState(i) # advance timeslider
    DrawPlots()
    SaveWindow()

# show time for last frame
#banner.text = str("Time = %5.2f days" % (time0-dt))
visit.View3DAttributes()
sys.exit()
