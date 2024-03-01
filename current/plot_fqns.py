from astropy.wcs import WCS
from astropy.io import fits
from astropy.nddata import Cutout2D
import astropy.units as un
from matplotlib.patches import Ellipse, Rectangle

def beam(hdu):
    """Creates a dictionary contating the information about the beam: b_major, b_minor, and b_pa. All units in
        degrees."""
    beam_dict = {};
    beam_dict = {"b_major": hdu.header['BMAJ'], "b_minor": hdu.header['BMIN'], "b_pa": hdu.header['BPA']};
    
    return(beam_dict)


def info(hdu, in_center_dict=None):
    """Creates a dictionary of information about the fits file: obs_code, obs_date, 
        obsra, obsdec, rms, min, and max."""
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
        raise TypeError("Type is" + str(hdu.__class___));
    
def load(filename):
    """Loads fits file into memory. Returns hdu, wcs, and hdu_data."""
    hdu = fits.open(filename)[0];
    hdu_data = hdu.data[0][0];
    wcs_tmp = WCS(hdu.header);
    wcs_tmp1 = wcs_tmp.dropaxis(3);
    wcs = wcs_tmp1.dropaxis(2);
    
    return(hdu_data, wcs, hdu);

def cut(data, in_wcs, hdu, center_dict, center_pos, box):
    """Cuts the data based on center position and size (box). 
        Returns data, wcs, hdu, and new center dict."""
    size = un.Quantity(box, un.pix);
    
    cutout = Cutout2D(data, position=center_pos, size=size, wcs=in_wcs);
    cutout_wcs = cutout.wcs.to_header();
    wcs_cut = WCS(cutout_wcs);
    cutout_data = cutout.data;
    
    center_tmp = cutout.center_cutout;
    #center_coord = wcs_cut.pixel_to_world(center_tmp[0],[1]);
    _ , new_center_dict = info(cutout_wcs, center_dict);
    new_center_dict |= {'centpxl_ra': center_tmp[0], 'centpxl_dec':center_tmp[1],
                      'center_ra':center_pos.ra.degree, 'center_dec':center_pos.dec.degree,
                       'npxl_ra': size[0].value, 'npxl_dec':size[1].value};
    
    return(cutout_data, wcs_cut, hdu, new_center_dict);

def scale(hdu, beam_dict, size):
    b_conv = (3600*1e3);
    beamsize = (beam_dict['b_major']*b_conv);
    sbar_scale = (size/beamsize);
    width = (beam_dict['b_major']*sbar_scale);

    return(width);

# ploting functions
def beam_plotter(beam_dict, center_dict, ax_name):
    """Returns a matplotlib Ellipse to plot on the image."""
    ra_offset = (center_dict['npxl_ra']/2); #35
    dec_offset = (center_dict['npxl_dec']/2 ); #2
    
    
    beam_ra = (center_dict['center_ra'] - ra_offset*center_dict['scale_ra'] -beam_dict['b_major']/2); 
    beam_dec = (center_dict['center_dec'] - dec_offset*center_dict['scale_dec'] + beam_dict['b_major']);
    
    beam = Ellipse((beam_ra, beam_dec), height=beam_dict['b_major'], width = beam_dict['b_minor'], 
                   angle=-beam_dict['b_pa'], transform=ax_name.get_transform('fk5'), edgecolor='w',
                   facecolor='w', alpha=0.85, lw=2); #angle needs to be neg as matplotlib angle is CCW
    return(beam)

def sbar_plotter(center_dict, beam_dict, sbar_width_ang, ax_name, finetune_ra=0, finetune_dec=0):
    """Returns scale bar object to be plotted, text, and the coordinates."""
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

def open(filename):
    """Loads the file into memory and creates the necessary dictionaries."""
    img_data, wcs, image = load(filename);
    
    img_info = info(image);
    beam_data = beam(image);
    
    allinone = (img_info | beam_data );
    return(img_data, wcs, image, allinone);
    
""" ef fit_plotter(name, ax_name):
    Returns the coordinates for the center of the fit, the fit shape as an Ellipse read from AIPS output.
    df = pd.read_csv(name, delimiter='\s\s+', engine='python', header=[0,1]);
    coords = SkyCoord(df.iat[0,3], df.iat[0,4], unit=(un.hourangle, un.degree), frame='fk5');
    theta = -df.iat[0,7];
    maj = ((df.iat[0,5])*un.mas).to(un.degree);
    min = ((df.iat[0,6])*un.mas).to(un.degree);
    fitra = coords.ra.degree;
    fitdec = coords.dec.degree;
    # fit = Ellipse((fitra, fitdec), height=maj, width=min, angle=theta,
    #               transform=ax_name.get_transform('fk5'), edgecolor='w', facecolor='w', 
    #               alpha=0.5, lw=2);
    fit = Ellipse((fitra, fitdec), height=maj.value, width = min.value, 
                   angle=theta, transform=ax_name.get_transform('fk5'), edgecolor='w',
                   facecolor='w', alpha=0.3, lw=2);
    return(coords, fit); """