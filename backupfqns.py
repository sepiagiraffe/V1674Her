#backup 2/18/25 before error analysis attempt
#-------------------------------------------------reading in files & whatnot fqns----------------------------------------------------------------------------
def beam(hdu):
    """input fits hdu"""
    beam_dict = {};
    beam_dict = {"b_major": hdu.header['BMAJ'], "b_minor": hdu.header['BMIN'], "b_pa": hdu.header['BPA']};
    #All units in degree
    
    return(beam_dict)


def info(hdu, in_center_dict=None):
    """inpt fits hdu"""
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
        print("The input is not a hdu object, it's:",str(hdu.__class__));

def load(filename):
    """input filename/path, will output hdu data, wcs, and hdu."""
    hdu = fits.open(filename)[0];
    hdu_data = hdu.data[0][0];
    wcs_tmp = WCS(hdu.header);
    wcs_tmp1 = wcs_tmp.dropaxis(3);
    wcs = wcs_tmp1.dropaxis(2);
    
    return(hdu_data, wcs, hdu)

def difmap_model_read(modelfit_file):
    """reads the output of difmap model files, path to file is needed."""
    header_names = ['Flux (Jy)', 'Flux (jy) Std', 'East (arcsec)', 'RA (arcsec) Std', 'North (arcsec)', 'Dec (arcsec) Std', 'Shape', 'R.A. (deg)', 'Dec (deg)',
        'Major FWHM (arcsec)', 'Major Std', 'Minor FWHM (arcsec)', 'Minor Std', 'Theta (deg)', 'Theta Std', 'Freq (Hz)', 'Spectral Index', 'Spec Indx Std' ];
    
    difmap_df = pd.read_csv(modelfit_file, header=None, sep=r'\s+', skiprows=[0,1,2], engine='python', names=header_names, skipinitialspace=True);
    # difmap_coords = 1;
    difmap_coords = SkyCoord(difmap_df.iloc[0,7], difmap_df.iloc[0,8], unit=(un.degree, un.degree),
                             obstime="J2000", frame='fk5');

    return(difmap_df, difmap_coords);

#-------------------------------------------------plotting & imaging fqns----------------------------------------------------------------------------
## Cutting function
def cut(data, in_wcs, center_dict, center_pos, box):
    """input data, wcs, center dictionary, center position and size of cutout."""
    size = un.Quantity(box, un.pix);
    
    cutout = Cutout2D(data, position=center_pos, size=size, wcs=in_wcs, copy=True);
    cutout_wcs = cutout.wcs.to_header();
    wcs_cut = WCS(cutout_wcs);
    cutout_data = cutout.data;
    
    center_tmp = cutout.center_cutout;
    # center_coord = wcs_cut.pixel_to_world(center_tmp[0],[1]);
    _ , new_center_dict = info(cutout_wcs, center_dict);
    # new_center_dict = {'centpxl_ra': center_coord.ra.degree, 'centpxl_dec':center_coord.dec.degree,
    #                   'center_ra': center_tmp[0], 'center_dec': center_tmp[1]};
    new_center_dict |= {'centpxl_ra': center_tmp[0], 'centpxl_dec':center_tmp[1],
                      'center_ra':center_pos.ra.degree, 'center_dec':center_pos.dec.degree,
                       'npxl_ra': size[0].value, 'npxl_dec':size[1].value};
    
    return(cutout_data, wcs_cut, new_center_dict);

##    # ploting functions
def scale(beam_dict, size):
    """input beam dictionary, and size of scale bar."""
    b_conv = (3600*1e3);
    beamsize = (beam_dict['b_major']*b_conv);
    sbar_scale = (size/beamsize);
    width = (beam_dict['b_major']*sbar_scale);

    return(width);

def beam_plotter(beam_dict, center_dict, ax_name, finetune_ra=0, finetune_dec=0):
    """finetune_ra (+) moves to the right, finetune_dec (+) moves down"""
    ra_offset = (center_dict['npxl_ra']/2 - 25  - finetune_ra); #- 25
    dec_offset = (center_dict['npxl_dec']/2 + finetune_dec )
    
    ra_offset_s = (ra_offset*center_dict['scale_ra']);
    dec_offset_s = (dec_offset*center_dict['scale_dec']);
   
    beam_ra = (center_dict['center_ra'] - ra_offset_s + beam_dict['b_major']/2); 
    beam_dec = (center_dict['center_dec'] - dec_offset_s + beam_dict['b_major']);
    
    beam = Ellipse((beam_ra, beam_dec), height=beam_dict['b_major'], width = beam_dict['b_minor'], 
                   angle=beam_dict['b_pa'], transform=ax_name.get_transform('fk5'), edgecolor='w',
                   facecolor='w', alpha=0.7, lw=2);
    return(beam);

def sbar_plotter(center_dict, beam_dict, sbar_width_ang, ax_name, finetune_ra=0, finetune_dec=0):
    """"Input center dictionary, beam dictionary, size of scale bar in mas image axis name
    additionally, you can fine tune the placement with fintune_ra & finetune_dec"""
    b_conv = (3600*1e3);
    beamsize = (beam_dict['b_major']*b_conv);
    sbar_scale = (sbar_width_ang/beamsize);
    
    sbarwidth = (beam_dict['b_major']*sbar_scale);
    
    sbar_ra_offset = (center_dict['npxl_ra']/2 - 20);
    sbar_dec_offset = (center_dict['npxl_dec']/2 - 15);
    
    sbar_rad = (center_dict['center_ra'] + sbar_ra_offset*center_dict['scale_ra']);
    sbar_dec = (center_dict['center_dec'] - sbar_dec_offset*center_dict['scale_dec']);
    
    
    sbartxt = str(sbar_width_ang) + ' mas';
  
    sbar_ra_txt = (sbar_ra_offset - 5*sbar_width_ang+ finetune_ra);
    sbar_dec_txt = (sbar_dec_offset- 7 + finetune_dec);
    
    sbar_txt_rad = (center_dict['center_ra'] +
                    sbar_ra_txt*center_dict['scale_ra']);
    sbar_txt_dec = (center_dict['center_dec'] - sbar_dec_txt*center_dict['scale_dec']);

    sbar = Rectangle((sbar_rad, sbar_dec), width = sbarwidth, height = beam_dict['b_minor']/10,
                    transform=ax_name.get_transform('fk5'), color='w', facecolor='w', lw=2);
                     
    return(sbar, sbartxt, sbar_txt_rad, sbar_txt_dec);


# -------------------------------------------------------------------general analysis functions ------------------------------------------------------------
def angular_size(distance: int, angular_size):
    """distance needs to be in pc and angular size in arcsec"""
    dist = distance*un.pc;
    input = angular_size*un.arcsec;
    theta = input.to(un.rad);
    Diam = np.tan(theta)*dist;
    Diam_AU = Diam.to(un.AU);  
    # print(f'{Diam:.3E}',f'{Diam_AU:.3E}');
    
    #some latex formating
    txt_di = '{D:.3E}';
    txt_di = txt_di.format(D=Diam.value);
    txt_dim = txt_di + "\," + "pc;";
    txt_da = '{F:.3E}';
    txt_da = txt_da.format(F=Diam_AU.value);
    txt_dam = txt_da + "\," + "AU";
    txt_diam = txt_dim + "\;\;\;\;" +txt_dam;
    display(Math(txt_diam));


    return  Diam, Diam_AU;


def angular_sizes(distance: int, angular_size, upper=0, lower=0):
    """distance needs to be in pc and angular size in arcsec"""
    if upper != 0:
        dist = distance*un.pc;
        dist_up = (distance+upper)*un.pc;
        dist_lo = (distance-lower)*un.pc;
        input = angular_size*un.arcsec;
        theta = input.to(un.rad);
        Diam_up = np.tan(theta)*dist_up;
        Diam_lo = np.tan(theta)*dist_lo;
        Diam = np.tan(theta)*dist;
        Diamlo_AU = Diam_lo.to(un.AU);  
        Diamup_AU = Diam_up.to(un.AU); 
        Diam_AU = Diam.to(un.AU); 
        print("Diameter (pc, AU):",f'{Diam:.3E}',f'{Diam_AU:.3E}',
              "Upper Diameter", f'{Diam_up:.3E}',f'{Diamup_AU:.3E}',
              "Lower Diameter", f'{Diam_lo:.3E}',f'{Diamlo_AU:.3E}' )
    elif upper ==0:    
        dist = distance*un.pc;
        input = angular_size*un.arcsec;
        theta = input.to(un.rad);
        Diam = np.tan(theta)*dist;
        Diam_AU = Diam.to(un.AU);  
        print(f'{Diam:.3E}',f'{Diam_AU:.3E}');
    
    return  Diam, Diam_AU;

def separation_mas(one, two, dist=0):
    if dist == 0:
        sep = one.separation(two);
        sep_mas = sep.to(un.mas);
        print(sep_mas);
    
    elif dist != 0:
        sep = one.separation(two);
        sep_mas = sep.to(un.mas);
        print(sep_mas);
        dist_pc = dist*un.pc;
        theta = sep_mas.to(un.rad);
        sep_pc = np.tan(theta)*dist_pc;
        sep_AU = sep_pc.to(un.AU);
        sep_cm = sep_pc.to(un.cm);
        print("Separation (pc):", f'{sep_pc:.3E}');
        print("Separation (AU):", f'{sep_AU:.3E}');
        print("Separation (cm):", f'{sep_cm:.3E}');

    return sep_mas

def magnetic_field(l, flux, Diameter, nu, K_o=40, f=1, i=0, si=-0.7):
    """path length (pc), flux (Jy), and Diameter (arcsec) 
    l must be an astropy unit object
    can specify K_o, f, i, or alpha (spectral index (si))
    assumes 10% of radius as the emitting region"""
    #Beck & Krause 2005

    l = (0.1*l.to(un.cm)).value;
    
    a = -1*si;
    S_nu = ((flux*un.Jy).cgs).value;
    D = (Diameter*un.arcsec).to(un.rad).value;
    
    I_nu = S_nu/D**2;
    E_p = ((938.257 *un.MeV).cgs).value;
    freq = ((nu*un.Hz).cgs).value;


    #constants 
    C3 = (1.86558e-23);
    C1 =(6.62428e18);
    # C2 = 0.25*C3* (a+5/3)/(a+1) * gamma(1+3a)*gamma((3*a+5)/6);
    C2 = 0.25*C3* (a+5/3)/(a+1) * gamma(1+3*a)*gamma((2*a+3)/4);
    # C2 = C2.decompose();
    C4 = np.cos(i)**(a+1);

    # exp = 1/(a+3);
    one = (4*np.pi*(2*a+1)*(K_o+1)*I_nu*(E_p**(1-2*a)));
    two = ((freq/(2*C1))**a);
    three = ((2*a-1)*C2*l*C4*f);
    B_eq = ((one*two)/three)**(1/(a+3));
    B_min = B_eq * ((a+1)/2)**(1/(a+3));
    # print(B_eq, "Gauss, K_0:", K_o, "path length (cm):", f'{l:.2E}');

    #some formatting in latex
    txt_B = "{B:.2E}"
    txt_B_var = "B_{eq}\: =";
    txt_B = txt_B.format(B=B_eq);
    txt_B = txt_B_var + txt_B  + "\,G";
    
    txt_l = "l\,=\,";
    txt_ln = "{L:.2E}";
    txt_ln = txt_ln.format(L=l);
    txt_l = txt_l + txt_ln + "\,cm"
    
    txt_K_0 = "K_0\,=\,";
    txt_K = '{K}';
    txt_K = txt_K.format(K=K_o);
    txt_K_0 = txt_K_0 + txt_K;
    
    txt_mega = txt_K_0 + ";\;\;\;" + txt_l + ";\;\;\;" + txt_B ;
    display(Math(txt_mega));

    return B_eq, B_min

def press_mag(B):
    """returns the magnetic field pressure in dyn/cm^2"""
    B_m = (B*un.G).to(un.T);
    P_b = (B_m)**2/(2*mu0);
    P_b = (P_b.value);

   
    #now to format in latex bc i hate myself
    #have to split the latex into 2 strings, otherwise the latex won't format
    txt0 = "P_B\,=\,";
    txt1 = "{x:.2E}";
    txt1 = txt1.format(x=P_b);
    txt2 = "\:\\frac{dyn}{cm^2}";
    txt = txt0 + txt1 + txt2;
    display(Math(txt));

    return P_b

def brightness_temp(theta_maj, theta_min, flux, freq):
    """ returns brightness temperature, needs theta in arcsec, flux in Jy, and freq"""
    rad_conv = (3600)**(-1)*(np.pi/180); #arcsec to rad
    jy = 10e-26;
    theta1 = theta_maj*rad_conv;
    theta2 = theta_min*rad_conv;
    k = 1.38e-23;
    c = 2.99e8;


    Tb = (c/freq)**2 * (flux*jy) * (2*k)**(-1) * ((4*np.log10(2))/(np.pi*theta1*theta2));

    print(f'{Tb:.2E}'+ " K");
    return Tb

def opacity_tb(T_b, T=1e4):
    """Returns the optical depth from brightness temperature"""
    tau = np.log(T_b/T -1);

    print(f'{tau:.2E}');
    return tau

def opacity_flux(F_trans, F_o):
    """returns the opacity from the ratio of transmitted and initial flux"""

    tau = np.log(F_trans/F_o);
    
    # print(f'{tau:.2E}');
    return tau

def emission_meas(tau, nu, T=10**4):
    """returns the emission measure in cm**6/pc"""
    A = tau/3.28e-7;
    B = (T/10e4)**1.35;
    C = (nu)**2.1;
    EM = A*B*C;

    # print(f'{EM:.2E}', "pc*cm**(-6)");
    return EM


def EM_FTR(F_trans, F_o, nu, l=0, T=10**4):
    """flux (Jy), l (pc); l is optional and will print n_e as a result"""
    tau = -1*opacity_flux(F_trans, F_o);
    EM = emission_meas(tau, nu, T);

    #latex shenanigans
    txt_t = "\\tau\,=\,";
    txt_tau = "{t:.2f}";
    txt_tau = txt_tau.format(t=tau);
    txt_op = txt_t + "\;" + txt_tau;
    
    # print("tau:",f'{tau:.2E}');
    # print("EM", f'{EM:.2E}')

    txt_EM = "EM\,=\,";
    txt_EM_v = "{x:.2E}";
    txt_EM_v = txt_EM_v.format(x=EM);
    txt_EM_un = "\\frac{pc}{cm^{6}}";
    txt_EM = txt_EM + txt_EM_v +"\,\,\," + txt_EM_un;

    if l != 0:
        
        #EM units: pc/cm^6
        l = (0.1*l).value;#untis: pc
        n_e = (EM/l)**(1/2); #units: 1/cm^3
        
        #latex
        txt_ne = "n_e\,=\,";
        txt_ne_v = "{x:.2E}";
        txt_ne_v = txt_ne_v.format(x=n_e);
        txt_ne_un = "\\frac{1}{cm^3}";
        txt_ne = txt_ne + txt_ne_v + "\;\;\;" + txt_ne_un;
        

        txt = txt_op + ";\;\;\;\;" + txt_EM + ";\;\;\;\;" + txt_ne;
        display(Math(txt));
    
    elif l ==0:
        tau = -1*opacity_flux(F_trans, F_o);
        EM = emission_meas(tau, nu, T);

        #latex
        txt = txt_op + ";\;\;\;\;" + txt_EM;
        display(Math(txt));

    return EM