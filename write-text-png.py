#!/usr/bin/env python

from PIL import Image, ImageDraw, ImageFont
import os,sys
import netCDF4
import numpy as np

def txt2img(text,file,FontSize=25):
    font_size = FontSize
    #fontPath = "/usr/share/fonts/dejavu-lgc/DejaVuLGCSansCondensed-Bold.ttf"
    fontPath = "/Library/Fonts/Verdana.ttf"
    sans16  =  ImageFont.truetype ( fontPath, FontSize )
    img = Image.open(file)
    draw = ImageDraw.Draw(img)                     # setup to draw on the main image
    draw.text((300,650), text,font=sans16, fill="white" )      # add some text to the main
    del draw
    img.save(file,"PNG",quality=100)

exp_name=sys.argv[1]
print 'Processing experiment',exp_name + '...'

# read non-dim time
data = np.loadtxt('time.txt')
time=data[:]
nf=len(time) # number of png
nf = 2
for fname in range(nf):
       print ' \n' + '==> ' + '  READING PNG FILES ...\n' + ''
       file = str(exp_name+'%04d.png' % (fname))
       s = str("txt2img('Time = %5.1f (days)',file)" % (time[fname]))
       eval(s)
       print fname

print ' \n' + '==> ' + '  DONE!\n' + ''
