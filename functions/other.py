#pylint: skip-file
#this is just a notes file
def opacity_tb(T_b, T=10**4):
    """Returns the optical depth from brightness temperature"""
    tau = np.log(T_b/T -1);

    print(f'{tau:.2E}');
    return tau


def opacity_flux(F_trans, F_o):
    """returns the opacity from the ratio of transmitted and initial flux"""

    tau = np.log(F_trans)/np.log(F_o);
    
    # print(f'{tau:.2E}');
    return tau

def EM_er(tau, dtau, nu, T=10**4, prints=True):

    EM = emission_meas(tau, nu, T);
    dEM = (1/3.28e-7)*(T/1e4)**1.35*nu**2.1*dtau;
    
    # print(EM, dEM);
    if prints is True:
        print(f'{EM:.2E}');
        print(f'{dEM:.2E}');

    return EM, dEM

def emission_meas(tau, nu, T=10**4):
    """returns the emission measure in cm**6/pc"""
    A = tau/3.28e-7;
    B = (T/1e4)**1.35;
    C = (nu)**2.1;
    EM = A*B*C;

    # print(f'{EM:.2E}', "pc*cm**(-6)");
    return EM

def n_e(F_trans, F_o, nu, l, T):

    tau = -1*opacity_flux(F_trans, F_o);
    EM = emission_meas(tau, nu, T);
    n_e = (EM/l)**(1/2); #units: 1/cm^3

    return tau, EM, n_e 

def opactity_flux_er(F_trans, F_o, dflux, prints=True):

    tau = (-1)*opacity_flux(F_trans, F_o); #this comes out negative so I'm just gonna multiply by -1
    dtau = np.log(F_trans)/(np.log(F_o)**2*F_o)*dflux;

    if prints is True:
        print(f'{tau:.2E}');
        print(f'{dtau:.2E}');


    return tau, dtau

def n_e_er(EM, dEM, l, dl, prints=True):

    n_e = (EM/l)**(1/2);

    
    dndEM = 1/2*(EM/l)**(-1/2)*l**(-1);
    dndl = -1/2*(EM/l)**(-1/2)*(EM/l**2);
    dn_e = np.sqrt((dndEM*dEM)**2+(dndl*dl)**2);

    if prints is True:
            
        print(f'{n_e:.2E}');
        print(f'{dn_e[0]:.2E}');
        print(f'{dn_e[1]:.2E}');

    return n_e, dn_e

def EM_FTR_er(F_trans, F_o, dflux, nu, l, dl, T=10**4, print_tex=False):

    tau, dtau = opactity_flux_er(F_trans, F_o, dflux, prints=False);

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


#emission stuff for free free absorption 

def EM(peak, rms, nu_difmap, Diam_df, beam, dflux, modelfit):
    
    l_in = Diam_df.iloc[0,3]; #the length of absorbing material is the minor axis size based on modelfit
    dll_in = Diam_df.iloc[0,4]; #lower error
    dlu_in = Diam_df.iloc[0,5]; #upper error

    beam_min = (beam['b_minor']*un.deg).to(un.arcsec).value;
    beam_maj = (beam['b_major']*un.deg).to(un.arcsec).value;
    #transmitted Flux (F_trans) calculation
    #3*rms*sqrt(size of the emission in # of beams)
    numb = beams(beam, modelfit);
    F_trans = 3*rms*np.sqrt(numb);

    conv = ((modelfit.iloc[0,9]*modelfit.iloc[0,11]*beam_min*beam_maj)*un.arcsec).to(un.mas).value;
    # F_trans = rms*3*conv; #we want 3 sigma detection
    nu = nu_difmap*1e-9;
    l = (l_in*.1).value; #10% of the minor axis
    dl_l = (dll_in*.1).value; 
    dl_u = (dlu_in*.1).value;
    dl = [dl_l, dl_u];
    
    op_list = EM_FTR_er(F_trans, peak, dflux, nu, l, dl);

    return op_list

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
        
        # a = (np.cos(theta)*np.sin(theta))**(-1)*(asize_er).value;
        a = (dist*(1+np.tan(theta)**2)*asize_er).value;
     

        for i in range(0,2):
            # b = (dist_er[i]/dist).value;
            
            b = (dist_er[i]*(np.tan(theta))).value;
            
            # print(np.sqrt(a**2+b**2)*Diam)
            er = (np.sqrt(a**2+b**2));
            
            # print(b, er)
            er_au = (er*un.pc).to(un.AU);
            
            error_pc.append(er);
            error_au.append(er_au);
    
        d_au = {'Major ang size': [Diam_AU[0]], 'majr lower er': [error_au[0][0]], 'major upper er': [error_au[1][0]], 
                'minor ang size': [Diam_AU[1]], 'minor lower er': [error_au[0][1]], 'minor upper er': [error_au[1][1]]};
        
        d_pc = {'Major ang size': [Diam[0]], 'majr lower er': [error_pc[0][0]], 'major upper er': [error_pc[1][0]], 
                'minor ang size': [Diam[1]], 'minor lower er': [error_pc[0][1]], 'minor upper er': [error_pc[1][1]]};
        
        d_au = pd.DataFrame(data=d_au);
        d_pc = pd.DataFrame(data=d_pc);

    return  d_au, d_pc

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
        
        # a = (np.cos(theta)*np.sin(theta))**(-1)*(asize_er).value;
        a = (dist*(1+np.tan(theta)**2)*asize_er).value;
     

        for i in range(0,2):
            # b = (dist_er[i]/dist).value;
            b = (dist_er[i]*(np.tan(theta))).value;
            # print(b)
            # print(np.sqrt(a**2+b**2)*Diam)
            er = (np.sqrt(a**2+b**2))*un.pc;
            
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
def press_mag_er(B, b_er_low, b_er_high, print_tex=False, tex=True):

    P_b = press_mag(B, tex=False);
    dP_bl = -(2*B/(np.pi*8))*b_er_low;
    dP_bu = (2*B/(np.pi*8))*b_er_high;
    

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