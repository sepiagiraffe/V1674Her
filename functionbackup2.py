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

def angsize_difmap(dist, difmap, dist_ler, dist_uer, print_tex=False):
    Diamb_au, Diamb_pc = angsize_er(dist, difmap.iloc[0,9], difmap.iloc[0,11], error=True, asize_maj_er=difmap.iloc[0,10],
                                    asize_min_er=difmap.iloc[0,12], dist_ler=dist_ler, dist_uer=dist_uer, pc=True, print_tex=print_tex);
    return Diamb_au, Diamb_pc

def angsize_er(distance: int, asize_maj, asize_min, error=False, asize_maj_er=0, asize_min_er=0, dist_ler=0, dist_uer=0, pc=False, print_tex=False):
    """distance needs to be in pc and angular size in arcsec"""

    if error is True:
        error_pc = [];
        error_au = [];

        dist = distance*un.pc;
        dist_er = [dist_ler, dist_uer]*un.pc;

        asize = [asize_maj, asize_min]*un.arcsec;
        theta = asize.to(un.rad);
        
        asize_er = [asize_maj_er, asize_min_er]*un.arcsec;
        asize_er = asize_er.to(un.rad);
        
        Diam = np.tan(theta)*dist;
        Diam_AU = Diam.to(un.AU); 
        
        a = (np.cos(theta)*np.sin(theta))**(-1)*(asize_er).value;
     

        for i in range(0,2):
        
            b = (dist_er[i]/dist).value;
            # print(b)
            # print(np.sqrt(a**2+b**2))
            er = (np.sqrt(a**2+b**2)*Diam);
            
            # print(b, er)
            er_au = er.to(un.AU);
            
            error_pc.append(er);
            error_au.append(er_au);
    
        #error format:
        #lower: maj, min ; upper: maj, min
        #error[0] gives lower: maj, min
        #error[0][0] gives lower, maj value & unit
        #now some latex formating for the print statement
        # error = error.to(un.AU);
        txt_majl_er = '{er:.2f}';
        txt_majl_er = txt_majl_er.format(er=-error_au[0][0].value);
        txt_minl_er = '{er:.2f}';
        txt_minl_er = txt_minl_er.format(er=-error_au[0][1].value);

        txt_maju_er = '{er:.2f}';
        txt_maju_er = txt_maju_er.format(er=error_au[1][0].value);
        txt_minu_er = '{er:.2f}';
        txt_minu_er = txt_minu_er.format(er=error_au[1][1].value);

        txt_maj = '{diam:.2f}';
        txt_maj = txt_maj.format(diam=Diam_AU[0].value);
        txt_min = '{diam:.2f}';
        txt_min = txt_min.format(diam=Diam_AU[1].value);
        
        display(Latex(f'${txt_maj}_{{{txt_majl_er}}}^{{{txt_maju_er}}}\; \mathrm{{AU}}$'));
        display(Latex(f'${txt_min}_{{{txt_minl_er}}}^{{{txt_minu_er}}}\; \mathrm{{AU}}$'));

        if pc is True:
            
            txt_majl_er = '{er:.2e}';
            txt_majl_er = txt_majl_er.format(er=-error_pc[0][0].value);
            txt_minl_er = '{er:.2e}';
            txt_minl_er = txt_minl_er.format(er=-error_pc[0][1].value);

            txt_maju_er = '{er:.2e}';
            txt_maju_er = txt_maju_er.format(er=error_pc[1][0].value);
            txt_minu_er = '{er:.2e}';
            txt_minu_er = txt_minu_er.format(er=error_pc[1][1].value);

            txt_maj = '{diam:.2e}';
            txt_maj = txt_maj.format(diam=Diam[0].value);
            txt_min = '{diam:.2e}';
            txt_min = txt_min.format(diam=Diam[1].value);
            
            display(Latex(f'${txt_maj}_{{{txt_majl_er}}}^{{{txt_maju_er}}}\; \mathrm{{pc}}$'));
            display(Latex(f'${txt_min}_{{{txt_minl_er}}}^{{{txt_minu_er}}}\; \mathrm{{pc}}$'));
        #give everything in AU
        
        d_au = {'Major ang size': [Diam_AU[0]], 'majr lower er': [error_au[0][0]], 'major upper er': [error_au[1][0]], 
                'minor ang size': [Diam_AU[1]], 'minor lower er': [error_au[0][1]], 'minor upper er': [error_au[1][1]]};
        
        d_pc = {'Major ang size': [Diam[0]], 'majr lower er': [error_pc[0][0]], 'major upper er': [error_pc[1][0]], 
                'minor ang size': [Diam[1]], 'minor lower er': [error_pc[0][1]], 'minor upper er': [error_pc[1][1]]};
        
        d_au = pd.DataFrame(data=d_au);
        d_pc = pd.DataFrame(data=d_pc);

        if print_tex:
            print(f'${txt_min}_{{{txt_minl_er}}}^{{{txt_minu_er}}}\; \mathrm{{AU}}$');

    
    elif error is False:
        dist = distance*un.pc;
        asize = [asize_maj, asize_min];
        input = asize*un.arcsec;
        theta = input.to(un.rad);
        Diam = np.tan(theta)*dist;
        Diam_AU = Diam.to(un.AU);  
        # print(f'{Diam:.3E}',f'{Diam_AU:.3E}');
        
        #some latex formating
        txt_di = '{D:.3E}';
        # print('tst')
        txt_di = txt_di.format(D=Diam[0].value);
        txt_dim = txt_di + "\," + "pc;";
        txt_da = '{F:.2f}';
        txt_da = txt_da.format(F=Diam_AU[0].value);
        txt_dam = txt_da + "\," + "AU";
        txt_diam = txt_dim + "\;\;\;\;" +txt_dam;

        
        txt_di_m = txt_di.format(D=Diam[1].value);
        txt_dim_m = txt_di_m + "\," + "pc;";
        txt_da_m = '{F:.2f}';
        txt_da_m = txt_da_m.format(F=Diam_AU[1].value);
        txt_dam_m = txt_da_m + "\," + "AU";
        txt_diam_m = txt_dim_m + "\;\;\;\;" +txt_dam_m;

        d_au = Diam_AU;
        d_pc = Diam;
        
        display(Latex(f'$\mathrm{{Major\;Axis=}}\; {txt_diam}$'));
        display(Latex(f'$\mathrm{{Minor\; Axis=}}\; {txt_diam_m}$'));

        if print_tex:
            print(f'Major = ${txt_diam}$');
            print(f'Minor = ${txt_diam_m}$')

    return  d_au, d_pc



def separation_mas(one, two, dist=0): #does this one need error bars
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


def B_field_er_difmap(difmap, Diam_pc, K_o=40):

    l = (Diam_pc.iloc[0,3]*.1);
    dl_in = ((Diam_pc.iloc[0,4]*.1).value, (Diam_pc.iloc[0,5]*.1).value);

    B_list = B_field_er(l, difmap.iloc[0,0], difmap.iloc[0,9], difmap.iloc[0,15], K_o=K_o,
                         error=True, dl=dl_in, dflux=difmap.iloc[0,1],dDiam=difmap.iloc[0,10]);
    
    return B_list


def B_field_er(l, flux, Diameter, nu, K_o=40, f=1, i=0, si=-0.7, error=False, dl=0, dflux=0, dDiam=0, print_tex=False):
    """path length (pc), flux (Jy), and Diameter (arcsec) 
    l must be an astropy unit object
    can specify K_o, f, i, or alpha (spectral index (si))
    assumes 10% of radius as the emitting region"""
    #Beck & Krause 2005
    if error is True:
        B_eq, B_min = magnetic_field(l, flux, Diameter, nu, K_o, f, i, si, tex=False);
        # print(B_eq, B_min);
        l = (l).to(un.cm).value;
    

        D = (Diameter*un.arcsec).to(un.rad).value;
        dflux = ((dflux*un.Jy).cgs).value;
        dDiam = (dDiam*un.arcsec).to(un.rad).value;
        dl = ((dl)*un.pc).to(un.cm).value;
    
        a = dflux/flux;
        b = 2*dDiam/D;
        c_l = dl[0]/l;
        c_u = dl[1]/l;
        # print(a)
       
        dB_eq_l = B_eq*np.sqrt(a**2+b**2+c_l**2);
        dB_eq_u = B_eq*np.sqrt(a**2+b**2+c_u**2);

        dB_min_l = B_min*np.sqrt(a**2+b**2+c_l**2);
        dB_min_u = B_min*np.sqrt(a**2+b**2+c_u**2);
         #some formatting in latex
        
        txt_l = "l=";
        txt_ln = "{L:.2E}";
        txt_ln = txt_ln.format(L=l);
        txt_l = txt_l + txt_ln + "\,cm"
        
        txt_K_0 = "K_0=";
        txt_K = '{K}';
        txt_K = txt_K.format(K=K_o);
        txt_K_0 = txt_K_0 + txt_K;
        
        txt_B = "{B:.3f}"
        txt_B = txt_B.format(B=B_eq);

        txt_B_min = "{B:.3f}"
        txt_B_min = txt_B_min.format(B=B_min);

        txt_l = '{l:.3f}';
        txt_l = txt_l.format(l=dB_eq_l);
        txt_u = '{u:.3f}';
        txt_u = txt_u.format(u=-dB_eq_u);

        txt_l_m = '{l:.3f}';
        txt_l_m = txt_l_m.format(l=dB_min_l);
        txt_u_m = '{u:.3f}';
        txt_u_m = txt_u_m.format(u=-dB_min_u);

        txt_mega = txt_K_0 + ";\;\;" ;

        # # B_list = [B_min, dB_min_l, dB_min_u, B_eq, dB_eq_l, dB_eq_u];
        # # cols = ['B min', 'B min er lower', 'B min er upper', 'B eq', 'B eq er lower', 'B eq er upper'];
        d = {'B min': [B_min], 'Bm lower er': [dB_min_l], 'Bm upper er': [dB_min_u], 'B eq': [B_eq],
             'B eq lower er': [dB_eq_l], 'B eq upper er': [dB_eq_u]};
        B_list = pd.DataFrame(data=d);

        display(Latex(f'${txt_mega}\; l={txt_ln}\;cm;\;B_{{eq}}={txt_B}_{{{txt_u}}}^{{{txt_l}}}\;(G);B_{{min}}={txt_B_min}_{{{txt_l_m}}}^{{{txt_u_m}}}\;(G)$'));

        if print_tex is True:
            
            print(f'${txt_mega}\; l={txt_ln}\;cm;\;B_{{eq}}={txt_B}_{{{txt_u}}}^{{{txt_l}}}\;(G);B_{{min}}={txt_B_min}_{{{txt_l_m}}}^{{{txt_u_m}}}\;(G)$');

    elif error is False:
        l = (l*un.pc).to(un.cm).value;
        B_eq, B_min = magnetic_field(l, flux, Diameter, nu, K_o, f, i, si, tex=True);

        d = {'B min': [B_min],  'B eq': [B_eq]};
        B_list = pd.DataFrame(data=d);

    return B_list

def magnetic_field(l, flux, Diameter, nu, K_o=40, f=1, i=0, si=-0.7, tex=True):
    """path length (pc), flux (Jy), and Diameter (arcsec) 
    l must be an astropy unit object
    can specify K_o, f, i, or alpha (spectral index (si))
    assumes 10% of radius as the emitting region"""
    #Beck & Krause 2005

    l = ((l).to(un.cm)).value;
    
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
    if tex is True:

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

def press_mage(b_list, tex=True, print_tex=False):
    """returns p_blist_min, p_blist_eq"""
    # press_mag_er(B, b_er_low, b_er_high, print_tex=False, tex=True)
    p_blist_min = press_mag_er(b_list.iloc[0][0], b_list.iloc[0][1], b_list.iloc[0][2], tex=tex, print_tex=print_tex);
    p_blist_eq = press_mag_er(b_list.iloc[0][3], b_list.iloc[0][4], b_list.iloc[0][5], tex=tex, print_tex=print_tex);

    return p_blist_min, p_blist_eq


def press_mag(B, tex=True, print_tex=False):
    """returns the magnetic field pressure in dyn/cm^2"""
    P_b = B**2/(8*np.pi);

    if tex is True:
        #now to format in latex bc i hate myself
        #have to split the latex into 2 strings, otherwise the latex won't format
        txt0 = "P_B\,=\,";
        txt1 = "{x:.2E}";
        txt1 = txt1.format(x=P_b);
        txt2 = "\:\\frac{dyn}{cm^2}";
        txt = txt0 + txt1 + txt2;
        display(Latex(f'${txt}$'));

    if print_tex is True:
        print(f'${txt}$');
    return P_b

def press_mag_er(B, b_er_low, b_er_high, print_tex=False, tex=True):

    P_b = press_mag(B, tex=False);
    dP_bl = -P_b*(2*b_er_low/B);
    dP_bu = P_b*(2*b_er_high/B);
    

    P_b_txt = '{P:.2e}';
    P_b_txt = P_b_txt.format(P=P_b);
    dP_bl_txt = f'{dP_bl:.2e}';
    dP_bu_txt = f'{dP_bu:.2e}';
    units = '\\frac{{dyn}}{{cm^2}}';
    display(Latex(f'${P_b_txt}_{{{dP_bl_txt}}}^{{{dP_bu_txt}}}\;{{{units}}}$'));

    if print_tex is True:
        print(f'${P_b_txt}_{{{dP_bl_txt}}}^{{{dP_bu_txt}}}\;{{{units}}}$');
    
    d = {'P_b': [P_b], 'P_b lower er': [dP_bl], 'P_B upper er': [dP_bu]};
    p_blist = pd.DataFrame(data=d);

    return p_blist

def dyn_press(n_e, v, prints=True):

    n_e_m = (n_e*(100)**3)*(1/un.m**3);
    rho = n_e_m*(m_p+m_e);

    P_dyn_SI = (1/2)*rho*v;
    P_dyn_cgs = P_dyn_SI.value/(1e-5*(100**2));
    
    if prints is True:
        print(f'{P_dyn_cgs:.2E}');


    return P_dyn_cgs

def dyn_press_er(n_e, v, dv, dn_e, print_tex=False):
    P_dyn = dyn_press(n_e, v, prints=False);


    a = (dn_e/n_e)**2;
    b = (2*dv/v)**2;
    dP = P_dyn*np.sqrt(a+b);
    dP[0] = -dP[0];

    P_dyn_txt = f'{P_dyn:.2E}';
    dPl_txt = f'{dP[0]:.2E}';
    dPu_txt = f'{dP[1]:.2E}';

    d = {'P (dyn/cm^2)': [P_dyn], 'P er (lower)': [dP[0]], 'P er (upper)': [dP[1]]};
    P_df = pd.DataFrame(data=d);

    p_dyn_txt = f'${P_dyn_txt}_{{{dPl_txt}}}^{{{dPu_txt}}}\\frac{{dyn}}{{cm^2}}$'
    display(Latex(p_dyn_txt));

    if print_tex is True:
        print(p_dyn_txt);

    return P_df

def dyn_press_list(op_list, v, dv, print_tex=False):
    P_df = dyn_press_er(op_list.iloc[0,4], v, dv, [op_list.iloc[0,5], op_list.iloc[0,6]], print_tex=print_tex);
     
    return P_df


def ideal_pres(n_e, T=10**4, prints=True):

    n_e = n_e*(1/un.cm**3);
    n_e_si = n_e.si;
    p = (n_e_si*k_B*(T*un.K));
    p_cgs = p.to(un.dyne/un.cm**2);
    
    if prints is True:
        print(f'{p_cgs:.2E}');

    return p_cgs

def ideal_pres_er(n_e, dn_e, T=10**4, print_tex=False):

    p_cgs = ideal_pres(n_e, T, prints=False);

    dp = p_cgs*(dn_e/n_e);
    dp[0] = -dp[0];

    d = {'P (dyne/cm^2)': [p_cgs], 'P er (lower)': [dp[0]], 'P er (upper)': [dp[1]]};
    p_df = pd.DataFrame(data=d);

    p_txt = f'{p_cgs.value:.2E}';
    dpl_txt = f'{dp[0].value:.2E}';
    dpu_txt = f'{dp[1].value:.2E}';

    tex = f'${p_txt}_{{{dpl_txt}}}^{{{dpu_txt}}}\;\\frac{{dyn}}{{cm^2}}$'
    display(Latex(tex));

    if print_tex is True:
        print(tex);
    
    return p_df

def ideal_pres_list(op_list, T=10**4, print_tex=False):
    p_df = ideal_pres_er(op_list.iloc[0,4], [op_list.iloc[0,5], op_list.iloc[0,6]], T, print_tex);

    return p_df

def brightness_temp_difmap(difmap, print_tex=False):

    Tb_df = brightness_temp_er(difmap.iloc[0,9], difmap.iloc[0,11], difmap.iloc[0,10], difmap.iloc[0,12], difmap.iloc[0,0],
                               difmap.iloc[0,1], difmap.iloc[0,15], print_tex=print_tex);
    return Tb_df


def brightness_temp(theta_maj, theta_min, flux, freq, tex=True):
    """ returns brightness temperature, needs theta in arcsec, flux in Jy, and freq"""
    rad_conv = (3600)**(-1)*(np.pi/180); #arcsec to rad
    jy = 10e-26;
    theta1 = theta_maj*rad_conv;
    theta2 = theta_min*rad_conv;
    k = 1.38e-23;
    c = 2.99e8;


    Tb = (c/freq)**2 * (flux*jy) * (2*k)**(-1) * ((4*np.log10(2))/(np.pi*theta1*theta2));
    Tb_txt = f'{Tb:.2e}';

    if tex is True:
       display(Latex(f'${Tb_txt}\;K$'));

    return Tb

def brightness_temp_er(theta_maj, theta_min, dtheta_maj, dtheta_min, flux, dflux, freq, print_tex=False):

    Tb = brightness_temp(theta_maj, theta_min, flux, freq, tex=False);

    dTb = Tb * np.sqrt((dflux/flux)**2+(dtheta_maj/theta_maj)**2+(dtheta_min/theta_min)**2);

    Tb_txt = f'{Tb:.2e}';
    dTb_txt = f'{dTb:.2e}';
    txt = f'${Tb_txt}\\pm\;{dTb_txt}\; K$';
    
    display(Latex(txt));

    d = {'Tb': [Tb], 'Tb er': [dTb]};
    Tb_df = pd.DataFrame(data=d);

    if print_tex is True:
        print(txt);

    return Tb_df

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

def EM_er(tau, dtau, nu, T=10**4, prints=True):

    EM = emission_meas(tau, nu, T);
    dEM = EM*(dtau/tau);
    
    # print(EM, dEM);
    if prints is True:
        print(f'{EM:.2E}');
        print(f'{dEM:.2E}');

    return EM, dEM

def emission_meas(tau, nu, T=10**4):
    """returns the emission measure in cm**6/pc"""
    A = tau/3.28e-7;
    B = (T/10e4)**1.35;
    C = (nu)**2.1;
    EM = A*B*C;

    # print(f'{EM:.2E}', "pc*cm**(-6)");
    return EM

def n_e(F_trans, F_o, nu, l, T):

    tau = -1*opacity_flux(F_trans, F_o);
    EM = emission_meas(tau, nu, T);
    n_e = (EM/l)**(1/2); #units: 1/cm^3

    return tau, EM, n_e 

def opactity_flux_er(F_trans, F_o, rms, prints=True):

    tau = (-1)*opacity_flux(F_trans, F_o); #this comes out negative so I'm just gonna multiply by -1
    dtau = (1/F_o)*rms;

    if prints is True:
        print(f'{tau:.2E}');
        print(f'{dtau:.2E}');


    return tau, dtau

def n_e_er(EM, dEM, l, dl, prints=True):

    n_e = (EM/l)**(1/2);

    a = (dEM/EM)**2;
    b = (dl/l)**2;
    sq = (1/2)*np.sqrt(a+b);

    dn_e = n_e*sq;

    if prints is True:
            
        print(f'{n_e:.2E}');
        print(f'{dn_e[0]:.2E}');
        print(f'{dn_e[1]:.2E}');

    return n_e, dn_e

def EM_FTR_er(F_trans, F_o, rms, nu, l, dl, T=10**4, print_tex=False):

    tau, dtau = opactity_flux_er(F_trans, F_o, rms, prints=False);

    tau_txt = f'{tau:.2}';
    dtau_txt = f'{dtau:.2}';

    EM, dEM = EM_er(tau, dtau, nu, T, prints=False);

    EM_txt = f'{EM:.2E}';
    dEM_txt = f'{dEM:.2E}';

    n_e, dn_e = n_e_er(EM, dEM, l, dl, prints=False);

    n_e_txt = f'{n_e:.2E}';
    dn_el_txt = f'{dn_e[0]:.2E}';
    dn_eu_txt = f'{dn_e[1]:.2E}';

    d = {'tau': [tau], 'dtau': [dtau], 'EM': [EM], 'dEM': [dEM], 'n_e': [n_e], 'dn_e lower': [dn_e[0]], 'dn_e upper': [dn_e[1]]};
    op_df = pd.DataFrame(data=d);
   
    tau_tex = f'$\\tau\; =\; {tau_txt}\pm{dtau_txt}\;$';
    EM_tex = (f'$EM\; =\; {EM_txt}\pm{dEM_txt}\;\\frac{{pc}}{{cm^6}}$');
    ne_tex = f'$n_e\;=\;{n_e_txt}^{{{dn_eu_txt}}}_{{{dn_el_txt}}}\;\\frac{{1}}{{cm^3}}$';

    display(Latex(tau_tex));
    display(Latex(EM_tex));
    display(Latex(ne_tex));

    if print_tex is True:
        print(tau_tex);
        print(EM_tex);
        print(ne_tex);

    return op_df

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

    d = {'tau': [tau], 'EM': [EM]}
    EM_df = pd.DataFrame(data=d);


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

        d_ne = {'tau': [tau], 'EM': [EM], 'n_e': [n_e]};
        EM_df = pd.DataFrame(data=d_ne);
    
    elif l ==0:
        # tau = -1*opacity_flux(F_trans, F_o);
        # EM = emission_meas(tau, nu, T);

        #latex
        txt = txt_op + ";\;\;\;\;" + txt_EM;
        display(Math(txt));

    return EM_df