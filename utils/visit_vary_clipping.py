import visit
import math
import csv

# initial position
c0 = visit.View3DAttributes()
c0.viewNormal = (0.904423, 0.0747664, 0.420034)
c0.focus = (625, 1000, -300.038)
c0.viewUp = (-0.409372, -0.125147, 0.903744)
c0.viewAngle = 30
c0.parallelScale = 1209.52
c0.nearPlane = -2430.67
c0.farPlane = 2430.67
c0.imagePan = (-0.0194483, -0.0190144)
c0.imageZoom = 1.45171
c0.perspective = 0
c0.eyeAngle = 2
c0.centerOfRotationSet = 0
c0.centerOfRotation = (625, 1000, -300.038)
c0.axis3DScaleFlag = 0
c0.axis3DScales = (1, 1, 1)
c0.shear = (0, 0, 1)

# set view to these values
visit.SetView3D(c0)

# number of snapshots
n=visit.TimeSliderGetNStates()

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
visit.SetOperatorOptions(ClipAtts, 2)

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
       SetActivePlots((2, 5))
       SetActivePlots(5)
       HideActivePlots()
    elif i<trho:
       ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1
       ClipAtts.plane1Status = 1
       ClipAtts.plane2Status = 0
       # update visit
       visit.SetOperatorOptions(ClipAtts, 2)

    visit.SetTimeSliderState(i) # advance timeslider
    visit.DrawPlots()
    visit.SaveWindow()

# show time for last frame
#banner.text = str("Time = %5.2f days" % (time0-dt))
visit.View3DAttributes()
