##Cutting function
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
