def angular_size(distance: int, angular_size, error=False, asize_maj_er=0, asize_min_er=0, dist_er=0):
    """distance needs to be in pc and angular size in arcsec"""

    if error:
        dist = distance*un.pc;
        dist_er = dist_er*un.pc;

        input = angular_size*un.arcsec;
        theta = input.to(un.rad);

        error = [asize_maj_er, asize_min_er];
        
       


        a = (np.cos(angular_size)*np.sin(angular_size))**(-1)*(asize_min_er);
        b = (dist_er/dist);

        a = (np.cos(angular_size)*np.sin(angular_size))**(-1)*(asize_er);
        b = (dist_er/dist);
        
        

        Diam = np.tan(theta)*dist;
        Diam_AU = Diam.to(un.AU); 

        lower_D_er = (np.sqrt(a**2+b**2)/Diam);
        upper_D_er = (np.sqrt(c**2+d**2)/Diam);
        

    elif error:
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


    return  Diam, Diam_AU