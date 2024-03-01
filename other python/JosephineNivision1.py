#!/usr/bin/env python
# coding: utf-8

# In[1]:


import astropy
from astropy.io import fits
import aplpy
import math
import matplotlib.pyplot as plt


# In[2]:


def drop2axes_AIPS(filename, outname):
    hdu = fits.open(filename)[0]
    for kw in "CTYPE", "CRVAL", "CRPIX", "CDELT", "CROTA", "NAXIS":
        for n in 3, 4:
            hdu.header.remove(f"{kw}{n}")
    fits.writeto(outname, hdu.data[0,0], hdu.header,overwrite=True)

def drop2axes_CASA(filename, outname):
    hdu = fits.open(filename)[0]
    for kw in "CTYPE", "CRVAL", "CRPIX", "CDELT", "CUNIT", "NAXIS":
        for n in 3, 4:
            hdu.header.remove(f"{kw}{n}")
    fits.writeto(outname, hdu.data[0,0], hdu.header,overwrite=True)
    


# Overall Image

# In[14]:


#open with atsropy first to input beam information & drop axes (this will write to a new file)
drop2axes_AIPS('V392PER.FITS','V392PER_D.FITS')
pic = fits.open('V392PER_D.FITS')
#make header easier to call
hdr = pic[0].header
#input beam shit, units is degrees (default in both AIPS and astropy)
hdr['BMAJ'] = 1.0044E-06 
hdr['BMIN'] = 5.2246E-07 
hdr['BPA'] = -1.52 


#if you want to check the header
print(repr(hdr))


# In[13]:


##now load with aplpy

fig = aplpy.FITSFigure(pic)

#cmap=plt.get_cmap('name_of_colormap')
#set max value (maxine) and rms 
maxine = 2.8590E-04
rms = 2.020E-05
#rms found in AIPS:
    #tvall, tvwin set box off source
    #imstat

#-----plot the image-----
fig.show_colorscale(vmin=rms,vmax=maxine) #vmin is frequency in our case. (aplpy default is velocity)
##,cmap=plt.get_cmap('CMRmap'), if you want to change the colormap to a matplotlib approved colormap
fig.add_colorbar()

#-----make the graph easier to read-----
fig.tick_labels.set_font(size='small', color='white')
fig.axis_labels.set_font(size='small', color='white')
#-----beam-------
fig.add_beam()
fig.beam.set_color('white')
#------recenter the image------
#center coords are from ds9 region
center_RA = 70.8390408
center_DEC = 47.3571771
width = 0.0000042
fig.recenter(center_RA, center_DEC, radius=width)
#some other fun stuff



#------contours-----
A = 3*rms
contours1 = [-1*A,0,1*A,A*math.sqrt(2),A*2,0.000121624,A*2*math.sqrt(2),0.00021376,A*4,
            .00025674,0.000275184,0.000295896,A*8,
            A*4*math.sqrt(2)*2]


fig.show_contour( levels= contours1,overlap=True)

#fig.save('/Users/mwiliam/Documents/Data/pics/V392PER.PNG')


# First Scan

# In[5]:


drop2axes_AIPS('V392PER_1.FITS','v392per1_1_D.fits')
pic1 = fits.open('v392per1_1_D.fits')
hdr1 = pic1[0].header
#input beam shit, units is degrees (default in both AIPS and astropy)
hdr1['BMAJ'] = 1.2051E-06 
hdr1['BMIN'] = 5.0349E-07
hdr1['BPA'] = 19.16 

print


# In[6]:


#now for the first half scans
fig1 = aplpy.FITSFigure(pic1)
#fix_aplpy_fits(fig1)
#get it to show some basic stuff, with easy to see colors
rms1 = 2.266E-05

fig1.show_colorscale(vmin=rms1,vmax=maxine)
#fig.add_grid()
fig1.tick_labels.set_font(size='small', color='white')
fig1.axis_labels.set_font(size='small', color='white')
#add beam
fig1.add_beam()
fig1.beam.set_color('white')
#recenter the image
#center coords are from ds9 region
center_RA = 70.8390411
center_DEC = 47.3571772

fig1.recenter(center_RA, center_DEC, radius=width)

A1 = 3 * rms1
contours2 = [-1*A1,0,1*A1,A1*math.sqrt(2),A1*2,A1*2*math.sqrt(2),0.00021376,A1*4,
            A1*8,
            A1*4*math.sqrt(2)*2]
fig1.show_contour(levels=contours2,layer = A)
fig1.add_colorbar()

#fig1.save('/Users/mwiliam/Documents/Data/pics/V392PER_1.PNG')


# Second Scan

# In[7]:


drop2axes_AIPS('V392PER_2.FITS','v392per1_2_D.fits')
pic2 = fits.open('v392per1_2_D.fits')
#make header easier to call
hdr2 = pic2[0].header
#input beam shit, units is degrees (default in both AIPS and astropy)
hdr2['BMAJ'] = 1.1986E-06 
hdr2['BMIN'] = 4.8704E-07
hdr2['BPA'] = -20.02 


# In[16]:


#now for the first half scans
fig2 = aplpy.FITSFigure(pic2)
#fix_aplpy_fits(fig2)
#get it to show some basic stuff, with easy to see colors
rms2 = 2.192E-05
A2 = 3 * rms2
contours3 = [-1*A2,0,1*A2,A2*math.sqrt(2),A2*2,A2*2*math.sqrt(2),0.00021376,A2*4,
            A2*8,A2*4*math.sqrt(2)*2]
cleo = 2.5449E-04
fig2.show_colorscale(vmin=rms2,vmax=cleo)

#vmax=2.5449E-04
#fig.add_grid()
fig2.tick_labels.set_font(size='small', color='white')
fig2.axis_labels.set_font(size='small', color='white')
#add beam
fig2.add_beam()
fig2.beam.set_color('white')
#recenter the image
#center coords are from ds9 region
center_RA = 70.8390411
center_DEC = 47.3571772
fig2.recenter(center_RA, center_DEC, radius=width)
#some other fun stuff

#fig2.show_contour(data= pic2,levels=contours,slices=[0,1])
#d = WCS('v392per1_1.fits')
#pic7 = WCS_fix(d)

#fig2.show_contour(levels=contours)
fig2.show_contour(data='v392per1_1_D.fits',levels=contours2)
fig2.add_colorbar()
print(maxine-2.5449E-04)


# In[17]:



fig2.save('/Users/mwiliam/Documents/Data/pics/V392PER_2_C_1.PNG')


# ### Same Beam, different times

# #### First Scan

# In[16]:


#open with atsropy first to input beam information & drop axes (this will write to a new file)
drop2axes_AIPS('v392per1b.fits','v392per1b_D.fits')
pic3 = fits.open('v392per1b_D.fits')
#make header easier to call
hdr3 = pic3[0].header
#input beam shit, units is degrees (default in both AIPS and astropy)
hdr3['BMAJ'] = 1.2051E-06
hdr3['BMIN'] = 1.2051E-06
hdr3['BPA'] = 0


#if you want to check the header
#print(repr(hdr3))


# In[59]:


fig3 = aplpy.FITSFigure(pic3)
fig3.show_colorscale()

fig3.tick_labels.set_font(size='small', color='white')
fig3.axis_labels.set_font(size='small', color='white')

fig3.add_beam()
fig3.beam.set_color('white')

center_RA = 70.8390411
center_DEC = 47.3571772
fig3.recenter(center_RA, center_DEC, radius=0.0000095)

rms3 = 9.625E-05
A3 = 3 * rms3
contours4 = [-1*A3,0,0.5*A3,0.75*A3,1*A3, 1.2*A3, 1.3*A3,2*A3,4*A3]
cleo3 = 3.8489E-04
fig3.show_colorscale(vmin=rms3,vmax=cleo3)

fig3.show_contour(levels=contours4)
#fig3.show_contour(returnlevels=True)
fig3.add_colorbar()
#fig3.save('/Users/mwiliam/Documents/Data/pics/V392PER1B.PNG')


# #### Second Scan 

# In[51]:


#open with atsropy first to input beam information & drop axes (this will write to a new file)
drop2axes_AIPS('v392per2b.fits','v392per2b_D.fits')
pic4 = fits.open('v392per2b_D.fits')
#make header easier to call
hdr4 = pic4[0].header
#input beam shit, units is degrees (default in both AIPS and astropy)
hdr4['BMAJ'] = 1.2051E-06
hdr4['BMIN'] = 1.2051E-06
hdr4['BPA'] = 0


# In[58]:


fig4 = aplpy.FITSFigure(pic4)
fig4.show_colorscale()

fig4.tick_labels.set_font(size='small', color='white')
fig4.axis_labels.set_font(size='small', color='white')

fig4.add_beam()
fig4.beam.set_color('white')
center_RA = 70.8390411
center_DEC = 47.3571772
fig4.recenter(center_RA, center_DEC, radius=0.0000095)

rms4 = 9.261E-05
A4 = 3 * rms4
contours5 = [-1*A4,0,0.5*A4,0.75*A4,1*A4, 1.2*A4, 1.3*A4,2*A4,4*A4]

cleo4 = 3.7548E-04
fig4.show_colorscale(vmin=rms4,vmax=cleo3)

fig4.show_contour(levels=contours5)
#fig4.show_contour(data='v392per1b_D.fits',levels=contours4)
fig4.add_colorbar()

fig4.save('/Users/mwiliam/Documents/Data/pics/V392PER2B.PNG')


# In[ ]:




