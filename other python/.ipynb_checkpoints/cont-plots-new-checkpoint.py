from astropy.wcs import WCS
import astropy
from astropy.io import fits
from astropy.coordinates import Angle
from astropy.utils.data import get_pkg_data_filename
# from regions import CirclePixelRegion, PixCoord
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator

import warnings
warnings.filterwarnings("ignore")

#%% load function
def load(filename):
    file = get_pkg_data_filename(filename)
    hdu = fits.open(file)[0]

    data = hdu.data
    wcs = WCS(hdu.header)
    
    return(data, wcs, hdu)


def find_center_and_scale(hdr):
    obsra, obsdec = hdr["OBSRA"], hdr["OBSDEC"]
    scale_ra, scale_dec = hdr["CDELT1"], hdr["CDELT2"] # deg/pixel
    # center_pixel_ra, center_pixel_dec = hdr["CRPIX1"], hdr["CRPIX2"]
    npxl_ra, npxl_dec = hdr["NAXIS1"], hdr["NAXIS2"]     
    return obsra, obsdec, scale_ra, scale_dec, npxl_ra, npxl_dec  # in deg

#%% cut function
def cut(data, wcs, hdu,center,box):
    size = u.Quantity(box, u.pix)
    cutout = Cutout2D(data,position=center,size=size,wcs=wcs)
    hdu.header.update(cutout.wcs.to_header())
    hdu.data = cutout.data

    wcs_cut = WCS(hdu.header)
    
    return(hdu.data,wcs_cut,hdu)


def cut2(data, wcs, hdu,center,box):
    size = u.Quantity(box, u.deg)
    cutout = Cutout2D(data,position=center,size=size,wcs=wcs)
    hdu.header.update(cutout.wcs.to_header())
    hdu.data = cutout.data

    wcs_cut = WCS(hdu.header)
    
    return(hdu.data,wcs_cut,hdu)

    
#%% find beam function
def find_beam(hdr):
    major = hdr["BMAJ"] # in degree
    minor = hdr["BMIN"]
    pa = hdr["BPA"] # position angle    
    return (major, minor, pa)


#%% -- 

file1 = '/home/tatush/Desktop/19035/fits_files/19035-Kband-Barray-full.fits' # K band
RoseroK_data , RoseroK_wcs , RoseroK_hdu = load(file1)
centerK_deg = SkyCoord(286.5067083333333,6.776722222222222, frame='fk5', unit='deg')
centerK_pix = RoseroK_wcs.world_to_pixel(centerK_deg)
# print(centerK_pix)
RoseroK_data , RoseroK_wcs , RoseroK_hdu = cut(RoseroK_data , RoseroK_wcs , 
                                               RoseroK_hdu, center=(1048,1224),box=(130,130))


file2 = '/home/tatush/Desktop/19035/fits_files/19035-Cband-Aarray-full.fits' # C band
RoseroC_data , RoseroC_wcs , RoseroC_hdu = load(file2)
centerC_deg = SkyCoord(286.5067083333333,6.776722222222222, frame='fk5', unit='deg')
centerC_pix = RoseroC_wcs.world_to_pixel(centerC_deg)
# print(centerC_pix)
RoseroC_data , RoseroC_wcs , RoseroC_hdu = cut(RoseroC_data , RoseroC_wcs , 
                                               RoseroC_hdu, (945,1116),box=(130,130))


mydata , mywcs , myhdu = load('/home/tatush/Desktop/19035/fits_files/19035-cont-0.013-im2000-pI-r0.5.pbcor.fits')
mydata = np.ma.masked_invalid(mydata)
rms=mydata.std()
center_deg = SkyCoord(286.5067083333333,6.776722222222222, frame='fk5', unit='deg')
center_pix = mywcs.world_to_pixel(center_deg)
# print(center_pix)
mydata , mywcs , myhdu = cut(mydata , mywcs , myhdu, center=(935,1020),box=(130,130))


## plot
fig = plt.figure(figsize=(15,15))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=mywcs)
    
    
    # -- axis label and ticks
ra = ax.coords[0]
dec = ax.coords[1]
ra.set_axislabel("RA (J2000)", minpad=0.8, fontsize=25)
dec.set_axislabel("Dec (J2000)", minpad=-1.0, fontsize=25)
ra.set_major_formatter('hh:mm:ss.ss')
ax.set_title('IRAS 19035+0641 A',fontsize=25)

    
ra.display_minor_ticks(True)
dec.display_minor_ticks(True)    
ax.tick_params(which='both',direction='in',color='white',length=10,width=2,
                   labelsize=20)
ax.tick_params(which='minor', length=5)

    # -- color plot
im=ax.imshow(mydata, vmin=-1e-9, vmax=0.00023, origin='lower', 
               cmap='inferno', transform=ax.get_transform(mywcs))
    

    # -- color bar
cbar = fig.colorbar(im,shrink=0.8)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(18)
cbar.ax.set_ylabel(r'$\mu$Jy/beam', fontsize=20)
cbar.ax.set_yticklabels(['0','50','100','150','200'])
    

plt.rcParams["lines.linewidth"] = 2
rmsK = 7e-6
levelsK = [rmsK*-3,rmsK*5,rmsK*10,rmsK*20,rmsK*30,rmsK*50,rmsK*57]      
ax.contour(RoseroK_data, levels=levelsK, colors='r', transform=ax.get_transform(RoseroK_wcs))
rmsC = 4e-6
levelsC = [rmsC*-5,rmsC*5,rmsC*10,rmsC*12,rmsC*20,rmsC*25]      
ax.contour(RoseroC_data, levels=levelsC, colors='cyan', transform=ax.get_transform(RoseroC_wcs))

# rms = 13.5e-6
# levels = [rms*-3,rms*3,rms*5,rms*10,rms*15]
# ax.contour(mydata, levels=levels, colors='r', transform=ax.get_transform(mywcs))


plt.rcParams["hatch.linewidth"] = 3
# Find image center and scale
obsra, obsdec, scale_ra, scale_dec, npxl_ra, npxl_dec = find_center_and_scale(myhdu.header)
major, minor, pa = find_beam(myhdu.header)
beam_ra  = obsra  - 0.85 *scale_ra*npxl_ra + major 
beam_dec = obsdec - 0.27 *scale_dec*npxl_dec + major

majorK, minorK, paK = find_beam(RoseroK_hdu.header)
beam2 = Ellipse((beam_ra, beam_dec), majorK, minorK, 90.0-paK, # BPA relative pos
                    transform=ax.get_transform('fk5'),
                    edgecolor='r', facecolor='none', hatch=r"//",lw=2)
ax.add_patch(beam2)

majorC, minorC, paC = find_beam(RoseroC_hdu.header)
# beam3 = Ellipse((beam_ra-3.7e-4, beam_dec-0.9e-5), majorC, minorC, 90.0-paC, # BPA relative pos
beam3 = Ellipse((beam_ra, beam_dec), majorC, minorC, 90.0-paC, # BPA relative pos
                    transform=ax.get_transform('fk5'),
                    edgecolor='c', facecolor='none', hatch="//",lw=3)
ax.add_patch(beam3)

# beam1 = Ellipse((beam_ra-7e-5, beam_dec), major, minor, 90.0-pa, # BPA relative pos
beam1 = Ellipse((beam_ra, beam_dec), major, minor, 90.0-pa, # BPA relative pos
                    transform=ax.get_transform('fk5'),
                    edgecolor='c', facecolor='w')
ax.add_patch(beam1)
    

## -- Scale bar
d = 2.3 # kpc
sb_length = (1e3) /(d * 1000) # I want a 1000 au scale bar
pix_scale = 0.013 # cell size
bsize = sb_length / pix_scale # size of bar in pix

xlim1, xlim2 = ax.get_xlim()
ylim1, ylim2 = ax.get_ylim()

aux = (xlim2-2-bsize+xlim2-2)/2.
ax.plot([xlim2-2-bsize,xlim2-2], 
                  [ylim2-3,ylim2-3], 
                  linewidth=3, c='w', alpha=1.0)

scale_bar_fontsize=16
scale_bar_text = '1000 au'

ax.text(aux, ylim2-8.5, scale_bar_text, 
                  fontsize=22,
                  horizontalalignment='center', 
                  color='w')

masers_ra = [286.5067794 , 286.50677574 ,286.50668674, 286.50662909,
 286.50658192 ,286.50652725 ] 
masers_dec = [ 6.77687944, 6.77684078 ,6.77671686 ,6.77665557 ,6.77664617,
 6.77652557]

# plt.plot(masers_ra,masers_dec,'w+', fillstyle='none', markersize=35, mew=2, transform=ax.get_transform('fk5'))

plt.text(286.50675,6.77688889,'A1',size=35,color='w', transform=ax.get_transform('fk5'))


# plt.savefig('19035-continuum.pdf',bbox_inches='tight',pad_inches=0.1)


#%%
## -- SECOND PLOT
fig = plt.figure(figsize=(15,15))
mydata , mywcs , myhdu = cut(mydata , mywcs , myhdu, center=(70,58),box=(70,70))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=mywcs)
    
    
    # -- axis label and ticks
ra = ax.coords[0]
dec = ax.coords[1]
ra.set_axislabel("RA (J2000)", minpad=0.8, fontsize=25)
dec.set_axislabel("Dec (J2000)", minpad=-1.0, fontsize=25)
ra.set_major_formatter('hh:mm:ss.ss')
ax.set_title('IRAS 19035+0641 A',fontsize=25)

    
ra.display_minor_ticks(True)
dec.display_minor_ticks(True)    
ax.tick_params(which='both',direction='in',color='white',length=10,width=2,
                    labelsize=20)
ax.tick_params(which='minor', length=5)


plt.rcParams["lines.linewidth"] = 1.5
    # -- color plot
im=ax.imshow(mydata, vmin=-1e-9, vmax=0.00023, origin='lower', 
                cmap='inferno', transform=ax.get_transform(mywcs))
    

    # -- color bar
cbar = fig.colorbar(im,shrink=0.8)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(18)
cbar.ax.set_ylabel(r'$\mu$Jy/beam', fontsize=20)
cbar.ax.set_yticklabels(['0','50','100','150','200'])

# Find image center and scale

beam_ra  = obsra  - 0.65 *scale_ra*npxl_ra + major 
beam_dec = obsdec - 0.18 *scale_dec*npxl_dec + major
beam = Ellipse((beam_ra, beam_dec), major, minor, 90.0-pa, # BPA relative pos
                    transform=ax.get_transform('fk5'),
                    edgecolor='w', facecolor='w',lw=3)
ax.add_patch(beam)
# print(mydata.max())

rms = 13.54e-6
# print(mydata.std())
levels = [rms*-2.5,rms*2.5,rms*5,rms*10,rms*15]
ax.contour(mydata, levels=levels, colors='w', transform=ax.get_transform(mywcs))

## -- Scale bar
d = 2.3 # kpc
sb_length = (500) /(d * 1000) # I want a 1000 au scale bar
pix_scale = 0.013 # cell size
bsize = sb_length / pix_scale # size of bar in pix

xlim1, xlim2 = ax.get_xlim()
ylim1, ylim2 = ax.get_ylim()

aux = (xlim2-2-bsize+xlim2-2)/2.
ax.plot([xlim2-2-bsize,xlim2-2], 
                  [ylim2-1.5,ylim2-1.5], 
                  linewidth=3, c='w', alpha=1.0)

scale_bar_fontsize=16
scale_bar_text = '500 au'

ax.text(aux, ylim2-4, scale_bar_text, 
                  fontsize=22,
                  horizontalalignment='center', 
                  color='w')

masers_ra = [286.50668674, 286.50662909,
 286.50658192] 
masers_dec = [ 6.77671686 ,6.77665557 ,6.77664617
 ]

# plt.plot(masers_ra,masers_dec,'w+', fillstyle='none', markersize=55, mew=2, transform=ax.get_transform('fk5'))


plt.text(286.50666667,6.77672222,'A2',size=35,color='w', transform=ax.get_transform('fk5'))
plt.text(286.50669583,6.77667222,'A3',size=35,color='w', transform=ax.get_transform('fk5'))

# plt.savefig('19035-continuum-zoom.pdf',bbox_inches='tight',pad_inches=0.1)
# print(rms)
        



