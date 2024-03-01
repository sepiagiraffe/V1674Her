#!/usr/bin/env python
# coding: utf-8

# In[1]:


import astropy
import aplpy


# In[2]:


#fix FITS images
def fix_aplpy_fits(aplpy_obj, dropaxis=2):
    #"""This removes the degenerated dimensions in APLpy 2.X...
    #The input must be the object returned by aplpy.FITSFigure().
    #`dropaxis` is the index where to start dropping the axis (by default it assumes the 3rd,4th place).
    #"""
    temp_wcs = aplpy_obj._wcs.dropaxis(dropaxis)
    temp_wcs = temp_wcs.dropaxis(dropaxis)
    aplpy_obj._wcs = temp_wcs

#To load fixed image:
#fig = aplpy.FITSFigure(pic)
#fix_aplpy_fits(fig)


# Overall Image

# In[3]:


#open with atsropy first to input beam information
pic = astropy.io.fits.open('V392PER.FITS')
#make header easier to call
hdr = pic[0].header
#input beam shit, units is degrees (default in both AIPS and astropy)
hdr['BMAJ'] = 1.0044E-06 
hdr['BMIN'] = 5.2246E-07 
hdr['BPA'] = -1.52 


#if you want to check the header
#print(repr(hdr))


# In[11]:


##now load with aplpy
fig = aplpy.FITSFigure(pic)
fix_aplpy_fits(fig)
#get it to show some basic stuff, with easy to see colors
fig.show_colorscale()
#fig.add_grid()
fig.tick_labels.set_font(size='small', color='white')
fig.axis_labels.set_font(size='small', color='white')
#add beam
fig.add_beam()
fig.beam.set_color('white')
#recenter the image
#center coords are from ds9 region
center_RA = 70.8390408
center_DEC = 47.3571771
fig.recenter(center_RA, center_DEC, radius=0.0000095)
#some other fun stuff
fig.add_colorbar()
#fig.colorbar.set_frame_color('white') doesn't work in this version :(
fig.show_colorscale(vmin=3.1062E-05,vmax=2.8590E-04) #vmin is frequency in our case. (aplpy default is velocity)
#so vmin = RMS noise, which was found in AIPS using IMEAN
# 2.9344E-05 is from hisogram
#3.1062E-05 is from pixels???
#3.502E-05 JY/BEAM is from all over????
#max value is 2.8590E-04
#min value is -1.0609E-04
C = 1.6704E-05 #this is a min
#contoursssss
#sqrt(2), 2, 2*sqrt(2), 4, 4*sqrt(2), 8
#contours = [ -(2)**(1/2)*C, -1*C, (2)**(1/2)*C, 2*C, 2*(2)**(1/2)*C, 4*C, 4*(2)**(1/2*C)]
contours = [-0.000106087,-4.07565e-05,2.4574e-05,8.99045e-05,0.000155235,0.000220565,0.000285896]
#print(contours)
fig.show_contour(levels=contours)
#print(levels)
fig.save('V392PER_E.PNG')


# First Scan

# In[5]:


pic1 = astropy.io.fits.open('v392per1_1.fits',verify='fix')
#make header easier to call
hdr = pic1[0].header
#input beam shit, units is degrees (default in both AIPS and astropy)
hdr['BMAJ'] = 1.2051E-06 
hdr['BMIN'] = 5.0349E-07
hdr['BPA'] = 19.16 


# In[10]:


#now for the first half scans
fig = aplpy.FITSFigure(pic1)
fix_aplpy_fits(fig)
#get it to show some basic stuff, with easy to see colors
fig.show_colorscale()
#fig.add_grid()
fig.tick_labels.set_font(size='small', color='white')
fig.axis_labels.set_font(size='small', color='white')
#add beam
fig.add_beam()
fig.beam.set_color('white')
#recenter the image
#center coords are from ds9 region
center_RA = 70.8390411
center_DEC = 47.3571772
#fig.recenter(center_RA, center_DEC, radius=0.0000100) #to get more of the top right emission
fig.recenter(center_RA, center_DEC, radius=0.0000095)
#some other fun stuff
#rms = 2.9300E-05
#max = 2.8590E-04
C = 1.3378E-04 #minimum 
contours = [-0.000133776,-7.09942e-05,-8.21233e-06,5.45695e-05,0.000117351,0.000180133,0.000242915]
fig.show_colorscale(vmin=2.9300E-05,vmax=2.8590E-04)
fig.show_contour(levels=contours)
fig.add_colorbar()
#fig.add_scalebar(length=.00000001,location='lower right')

# ------------ Spots ----------------
S_TL_RA = 70.8390373
S_BL_RA = 70.8390448
S_TL_DEC = 47.3571824
S_BL_DEC = 47.3571713
S_BL_R = 0.0000028
S_TL_R = 0.0000032


# Second Scan

# In[7]:


pic2 = astropy.io.fits.open('v392per1_2.fits',verify='fix')
#make header easier to call
hdr = pic2[0].header
#input beam shit, units is degrees (default in both AIPS and astropy)
hdr['BMAJ'] = 1.1986E-06 
hdr['BMIN'] = 4.8704E-07
hdr['BPA'] = -20.02 


# In[12]:


#now for the first half scans
fig = aplpy.FITSFigure(pic2)
fix_aplpy_fits(fig)
#get it to show some basic stuff, with easy to see colors
fig.show_colorscale()
#fig.add_grid()
fig.tick_labels.set_font(size='small', color='white')
fig.axis_labels.set_font(size='small', color='white')
#add beam
fig.add_beam()
fig.beam.set_color('white')
#recenter the image
#center coords are from ds9 region
center_RA = 70.8390411
center_DEC = 47.3571772
fig.recenter(center_RA, center_DEC, radius=0.0000095)
#some other fun stuff
C = 1.3297E-04 #minimum 
contours = [-0.000133776,-7.09942e-05,-8.21233e-06,5.45695e-05,0.000117351,0.000180133,0.000242915]
fig.show_colorscale(vmin=3.0551E-05,vmax=2.5449E-04)
fig.show_contour()
fig.add_colorbar()

#------Spots--------
S_TL_RA_2 = 70.8390448
S_BL_RA_2 = 70.8390498
S_TL_DEC_2 = 47.3571713
S_BL_DEC_2 = 47.3571740
S_BL_R_2 = 0.0000030
S_TL_R_2 = 0.0000028
fig.show_circles(S_TL_RA,S_TL_DEC,S_TL_R,color='white')
fig.show_circles(S_BL_RA,S_BL_DEC,S_BL_R,color='white')


# In[ ]:




