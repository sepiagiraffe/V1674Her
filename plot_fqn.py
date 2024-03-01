import astropy
from astropy.wcs import WCS
from astropy import wcs
from astropy.io import fits, ascii
from astropy.coordinates import Angle, SkyCoord
from astropy.nddata import Cutout2D
import astropy.units as un 


import pandas as pd


from matplotlib.patches import Ellipse, Rectangle
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator, FixedFormatter

def beam(hdu):
    beam_dict = {};
    beam_dict = {"b_major": hdu.header['BMAJ'], "b_minor": hdu.header['BMIN'], "b_pa": hdu.header['BPA']};
    #All units in degree
    
    return(beam_dict)


def info(hdu, in_center_dict=None):
    
    if (str(hdu.__class__) == "<class 'astropy.io.fits.hdu.image.PrimaryHDU'>"):
        
        info_dict = {};
        info_dict = { "obs_code": hdu.header['OBSERVER'], "obs_date": hdu.header['DATE-OBS'],
                 "obsra": hdu.header['OBSRA'], "obsdec": hdu.header['OBSDEC'],
                 "rms": hdu.header['NOISE'], "min": hdu.header['DATAMIN'], "max": hdu.header['DATAMAX']};
        center_dict = {};
        center_dict = {"centpxl_ra": hdu.header['CRPIX1'], "centpxl_dec": hdu.header['CRPIX2'],
                "center_ra": hdu.header['CRVAL1'], "center_dec": hdu.header['CRVAL2'],
                "scale_ra": hdu.header['CDELT1'], "scale_dec": hdu.header['CDELT2'],
                "npxl_ra": hdu.header['NAXIS1'], "npxl_dec": hdu.header['NAXIS2']};
        return(info_dict, center_dict);
        
    elif (str(hdu.__class__) == "<class 'astropy.io.fits.header.Header'>"):
        
        info_dict = {};
        info_dict = { "obs_date": hdu['DATE-OBS']};
        center_dict = in_center_dict;
        center_dict |= {"centpxl_ra": hdu['CRPIX1'], "centpxl_dec": hdu['CRPIX2'],
                "center_ra": hdu['CRVAL1'], "center_dec": hdu['CRVAL2'],
                "scale_ra": hdu['CDELT1'], "scale_dec": hdu['CDELT2']};
        return(info_dict, center_dict);
    else:
        print("Not work!",str(hdu.__class___));

def load(filename):
    hdu = fits.open(filename)[0];
    hdu_data = hdu.data[0][0];
    wcs_tmp = WCS(hdu.header);
    wcs_tmp1 = wcs_tmp.dropaxis(3);
    wcs = wcs_tmp1.dropaxis(2);
    
    return(hdu_data, wcs, hdu)

def open(filename):
    img_data, wcs, image = load(filename);
    
    img_info, center_data = info(image);
    beam_data = beam(image);
    
    allinone = (img_info | beam_data | center_data);
    return(img_data, wcs, image, allinone)

## Cutting function
def cut(data, in_wcs, hdu, center_dict, center_pos, box):
    
    size = un.Quantity(box, un.pix);
    
    cutout = Cutout2D(data, position=center_pos, size=size, wcs=in_wcs, copy=True);
    cutout_wcs = cutout.wcs.to_header();
    wcs_cut = WCS(cutout_wcs);
    cutout_data = cutout.data;
    
    center_tmp = cutout.center_cutout;
    center_coord = wcs_cut.pixel_to_world(center_tmp[0],[1]);
    _ , new_center_dict = info(cutout_wcs, center_dict);
    # new_center_dict = {'centpxl_ra': center_coord.ra.degree, 'centpxl_dec':center_coord.dec.degree,
    #                   'center_ra': center_tmp[0], 'center_dec': center_tmp[1]};
    new_center_dict |= {'centpxl_ra': center_tmp[0], 'centpxl_dec':center_tmp[1],
                      'center_ra':center_pos.ra.degree, 'center_dec':center_pos.dec.degree,
                       'npxl_ra': size[0].value, 'npxl_dec':size[1].value};
    
    return(cutout_data, wcs_cut, hdu, new_center_dict);

##    # ploting functions
def scale(hdu, beam_dict, size):
    b_conv = (3600*1e3);
    beamsize = (beam_dict['b_major']*b_conv);
    sbar_scale = (size/beamsize);
    width = (beam_dict['b_major']*sbar_scale);

    return(width);

def beam_plotter(beam_dict, center_dict, ax_name):
    
    ra_offset = (center_dict['npxl_ra']/2 + 10);
    dec_offset = (center_dict['npxl_dec']/2 - 5);
    
    beam_ra = (center_dict['center_ra'] - ra_offset*center_dict['scale_ra'] + beam_dict['b_major']); 
    beam_dec = (center_dict['center_dec'] - dec_offset*center_dict['scale_dec'] + beam_dict['b_major']);
    
    beam = Ellipse((beam_ra, beam_dec), height=beam_dict['b_major'], width = beam_dict['b_minor'], 
                   angle=beam_dict['b_pa'], transform=ax_name.get_transform('fk5'), edgecolor='w',
                   facecolor='w', alpha=0.7, lw=2);
    return(beam, beam_ra, beam_dec);

def sbar_plotter(center_dict, beam_dict, sbar_width_ang, ax_name, finetune_ra=0, finetune_dec=0):
    
    b_conv = (3600*1e3);
    beamsize = (beam_dict['b_major']*b_conv);
    sbar_scale = (sbar_width_ang/beamsize);
    
    sbarwidth = (beam_dict['b_major']*sbar_scale);
    
    sbar_ra_offset = (center_dict['npxl_ra']/2 - 20);
    sbar_dec_offset = (center_dict['npxl_dec']/2 - 15);
    
    sbar_rad = (center_dict['center_ra'] + sbar_ra_offset*center_dict['scale_ra']);
    sbar_dec = (center_dict['center_dec'] - sbar_dec_offset*center_dict['scale_dec']);
    
    
    sbartxt = str(sbar_width_ang) + ' mas';
  
    sbar_ra_txt = (sbar_ra_offset - 5.25*sbar_width_ang+ finetune_ra);
    sbar_dec_txt = (sbar_dec_offset- 7 + finetune_dec);
    
    sbar_txt_rad = (center_dict['center_ra'] +
                    sbar_ra_txt*center_dict['scale_ra']);
    sbar_txt_dec = (center_dict['center_dec'] - sbar_dec_txt*center_dict['scale_dec']);

    sbar = Rectangle((sbar_rad, sbar_dec), width = sbarwidth, height = beam_dict['b_minor']/10,
                    transform=ax_name.get_transform('fk5'), color='w', facecolor='w', lw=2);
                     
    return(sbar, sbartxt, sbar_txt_rad, sbar_txt_dec);

def difmap_model_read(modelfit_file):

    header_names = ['Flux (Jy)', 'Flux (jy) Std', 'East (arcsec)', 'RA (arcsec) Std', 'North (arcsec)', 'Dec (arcsec) Std', 'Shape', 'R.A. (deg)', 'Dec (deg)',
        'Major FWHM (arcsec)', 'Major Std', 'Minor FWHM (arcsec)', 'Minor Std', 'Theta (deg)', 'Theta Std', 'Freq (Hz)', 'Spectral Index', 'Spec Indx Std' ];
    
    difmap_df = pd.read_csv(modelfit_file, header=None, sep=r'\s+', skiprows=[0,1,2], engine='python', names=header_names, skipinitialspace=True);
    # difmap_coords = 1;
    difmap_coords = SkyCoord(difmap_df.iloc[0,7], difmap_df.iloc[0,8], unit=(un.degree, un.degree),
                             obstime="J2000", frame='fk5');

    return(difmap_df, difmap_coords);
