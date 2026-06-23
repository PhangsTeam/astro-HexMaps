import numpy as np
import pandas as pd
import os.path
# from scipy import stats
from astropy import units as au
from astropy.stats import median_absolute_deviation, mad_std
from astropy.table import Table, Column


from structure_addition import *
from shuffle_spec import *
from mom_computer import get_mom_maps

def construct_mask(ref_line, this_data, SN_processing):
    """
    Function to construct the mask based on high and low SN cut
    """
    ref_line_data = this_data["SPEC_"+ref_line]
    n_pts = np.shape(ref_line_data)[0]
    n_chan = np.shape(ref_line_data)[1]

    line_vaxis = this_data.meta['SPEC_VCHAN0']+(np.arange(n_chan)-(this_data.meta['SPEC_CRPIX']-1))*this_data.meta['SPEC_DELTAV']

    line_vaxis = line_vaxis.to(au.km/au.s) #to km/s
    #Estimate rms
    rms = median_absolute_deviation(ref_line_data, axis = None, ignore_nan = True)
    rms = median_absolute_deviation(ref_line_data[np.where(ref_line_data<3*rms)], ignore_nan = True)

    #First create a rough mask
    mask_rough = (ref_line_data) < 3 * rms
    masked_cube = np.where(mask_rough, ref_line_data, np.nan)
    med_mask = np.nanmedian(masked_cube, axis=1)

    # Median absolute deviation along z (ignoring NaNs)
    mad_mask = np.nanmedian(np.abs(masked_cube - med_mask[:,None]), axis=1)

    # Mask each spectrum
    low_tresh, high_tresh = SN_processing[0]* mad_mask[:, None], SN_processing[1]* mad_mask[:, None]    
    mask = np.array(ref_line_data > high_tresh , dtype = int)
    low_mask = np.array(ref_line_data > low_tresh , dtype = int)

    mask = mask & (np.roll(mask, 1,1) | np.roll(mask,-1,1))

    #remove spikes along spectral axis:
    mask = np.array((mask + np.roll(mask, 1, 1) + np.roll(mask, -1, 1))>=3, dtype = int)
    low_mask = np.array((low_mask + np.roll(low_mask, 1, 1) + np.roll(low_mask, -1, 1))>=3, dtype = int)

    #remove spikes along spatial axis:
    #mask = np.array((mask + np.roll(mask, 1, 0) + np.roll(mask, -1, 0))>=3, dtype = int)
    #low_mask = np.array((low_mask + np.roll(mask, 1, 0) + np.roll(low_mask, -1, 0))>=3, dtype = int)

    #expand to cover all > 2sigma that have a 2-at-4sigma core
    for kk in range(5):
        mask = np.array(((mask + np.roll(mask, 1, 1) + np.roll(mask, -1, 1)) >= 1), dtype = int)*low_mask

    #expand to cover part of edge of the emission line
    for kk in range(2):
        mask = np.array(((mask + np.roll(mask, 1, 1) + np.roll(mask, -1, 1)) >= 1), dtype = int)
    
    # Derive the ref line mean velocity
    line_vmean = np.zeros(n_pts)*np.nan * au.km / au.s
    mask *=au.dimensionless_unscaled
    for jj in range(n_pts):
        line_vmean[jj] = np.nansum(line_vaxis * ref_line_data[jj,:]*mask[jj,:])/ \
                       np.nansum(ref_line_data[jj,:]*mask[jj,:])

    return mask, line_vmean, line_vaxis

def dist(ra, dec, ra_c, dec_c):
    return np.sqrt((ra-ra_c)**2+(dec-dec_c)**2)
    
def process_spectra(source_list,
                    lines_data,
                    fname,
                    shuff_axis,
                    run_success,
                    ref_line_method,
                    SN_processing,
                    strict_mask,
                    input_mask = None,
                    use_input_mask = False,
                    hfs_data = None,
                    use_hfs_lines = False,
                    mom_calc = [3, 3, "fwhm"],
                    just_source = None,
                    use_noise_vel_ranges = False,
                    ):
    """
    :param lines_data:        Pandas DataFrame which is the cubes_list.txt.
    :param use_noise_vel_ranges: If True, read the SPEC_NOISE_MASK column from the
                                 database and pass it to get_mom_maps so that the RMS
                                 is estimated only over the explicitly defined noise
                                 velocity window(s) rather than all off-signal channels.
    """
    
    n_sources = len(source_list)
    n_lines = len(lines_data["line_name"])
    if ref_line_method in list(lines_data["line_name"]):
        #user defined reference line
        ref_line = ref_line_method.upper()
    else:
        ref_line = lines_data["line_name"][0].upper()

    for ii in range(n_sources):

        #if the run was not succefull, don't do processing of the data
        if not run_success[ii]:
            continue

        this_source = source_list[ii]
        if not just_source is None:
            if just_source != this_source:
                continue

        # print("----------------------------------")
        # print(f'Source: {this_source}')
        # print("----------------------------------")

        this_data = Table.read(fname[ii])
        tags = this_data.keys()
        n_chan = np.shape(this_data["SPEC_"+ref_line])[1]

        #--------------------------------------------------------------
        #  Build a mask based on reference line(s)
        #--------------------------------------------------------------
        
        if use_input_mask:
            # print(f'{"[INFO]":<10}', 'Use input mask.')
            # check that input mask is provided
            if len(input_mask) == 0:
                print(f'{"[ERROR]":<10}', f'No mask provided!')

            # use input mask
            # mask = this_data["SPEC_VAL_MASK"]
            mask = this_data[f'SPEC_{input_mask["mask_name"][0].upper()}']

            # clean up database
            del this_data[f'SPEC_{input_mask["mask_name"][0].upper()}']

            # take reference velocity and vaxis from reference line
            _, ref_line_vmean, ref_line_vaxis = construct_mask(ref_line, this_data, SN_processing)

        else:
            print(f'{"[INFO]":<10}', 'Build mask from prior line(s).')
            # Use function for mask
            mask, ref_line_vmean, ref_line_vaxis = construct_mask(ref_line, this_data, SN_processing)
            this_data["SPEC_MASK_"+ref_line]= Column(mask, unit=au.dimensionless_unscaled, description=f'Velocity-integration mask for {ref_line}')

            # check which lines are used as priors
            # n_mask = 0
            line_names = [str(l) for l in lines_data['line_name']]
            if ref_line_method in ["first"]:
                n_mask = 0
                print(f'{"[INFO]":<10}', f'Using first line as prior: {line_names[0]}.')
            elif ref_line_method in ["all"]:
                n_mask = n_lines
                print(f'{"[INFO]":<10}', f'All lines used as prior: {line_names}.')
            elif isinstance(ref_line_method, int):
                n_mask = np.min([n_lines,ref_line_method])
                print(f'{"[INFO]":<10}', f'Using first {n_mask+1} lines as prior: {line_names[:n_mask+1]}.')

            # combine masks of all priors if more than one line is used
            if n_mask>0:
                # for n_mask_i in range(1,n_mask): 
                for n_mask_i in range(1,n_mask+1):  
                    line_i = lines_data["line_name"][n_mask_i].upper()
                    mask_i, _, _ = construct_mask(line_i, this_data, SN_processing)
                    this_data["SPEC_MASK_"+line_i]= Column(mask_i, unit=au.dimensionless_unscaled, description=f'Velocity-integration mask for {line_i}')

                    # add mask to existing mask
                    mask = mask.value.astype(int) | mask_i.value.astype(int)
                    mask *=au.dimensionless_unscaled

            # combine with HI mask
            if ref_line_method in ["ref+HI"]:
                if "hi" not in list(lines_data["line_name"]):
                    print(f'{"[WARNING]":<10}', 'HI not in PyStructure. Skipping.')
                else:
                    mask_hi, ref_line_vmean_hi, ref_line_vaxis_hi = construct_mask("HI", this_data, SN_processing)
                    
                    # add mask to existing mask
                    mask = mask.value.astype(int) | mask_hi.value.astype(int)
                    mask *=au.dimensionless_unscaled
                    
                    rgal = this_data["rgal_r25"]
                    n_pts = len(this_data["rgal_r25"])
                    vmean_comb = np.zeros(n_pts)*np.nan
                    for jj in range(n_pts):
                        if rgal[jj]<0.23:
                            vmean_comb[jj] = ref_line_vmean[jj]
                        else:
                            vmean_comb[jj] = ref_line_vmean_hi[jj]
                    ref_line_vmean = vmean_comb

            if strict_mask:
                """
                Make sure that spatially we do not have only connected pixels
                """
                ra, dec = this_data["ra_deg"], this_data["dec_deg"]
                for jj in range(n_chan):
                    mask_spec = mask[:,jj]
                    mask_labels=np.zeros_like(mask_spec)
                    sep = this_data["beam_as"]/3600/2
                    label=1
                    for n in range(len(mask_labels)):
                        if mask_labels[n]==0:
                            if mask_spec[n]==0:
                                mask_labels[n]=-99
                                continue
            
                            dist_array=dist(ra, dec, ra[n], dec[n])
                            #check out neighbours
                            idx_neigh=np.where(abs(dist_array-sep)<0.1*this_data.meta["beam_as"].to(au.deg))
                            #check if labels have already been given (except 0 or -99)
                            labels_given=np.unique(mask_labels[idx_neigh])
                            index = labels_given[labels_given>0]
                            if len(index)>0:
                                mask_labels[n]=index[0]
                                if len(index)>1:
                                    for i in range(len(index)-1):
                                        mask_labels[mask_labels==index[i+1]]=index[0]
                            else:
                                mask_labels[n]=label
                                label+=1
                    labels = np.unique(mask_labels)
                    for lab in labels:
                        if lab <0:
                            continue
                        if len(mask[:,jj][np.where(mask_labels==lab)])<5:
                            mask[:,jj][np.where(mask_labels==lab)]=0

        #-------------------------------------------------------------------
        # OPTIONAL: Modified mask for hyperfine structure lines
        #-------------------------------------------------------------------

        if use_hfs_lines:

            # check that input hyperfine structure data is provided
            if len(hfs_data) == 0:
                print(f'{"[ERROR]":<10}', 'No hyperfine structure file provided!')

            # get set of hfs lines
            lines_hfs = list(set(hfs_data['hfs_name']))

            # get associated reference rest frequency
            line_restfreq_list = []
            for line in lines_hfs:
                # select columns of given line
                idx_cols = hfs_data['hfs_name'] == line
                # get list of frequencies for given line
                restfreqs = [f * au.Unit(str(u)) for f, u in zip(hfs_data['hfs_ref_freq'][idx_cols], hfs_data['unit'][idx_cols])]
                # append to list
                line_restfreq_list.append(restfreqs)

            # get set of hyperfine structure frequencies per line
            hfs_freq_list = []
            for line in lines_hfs:
                # select columns of given line
                idx_cols = hfs_data['hfs_name'] == line
                # get list of frequencies for given line
                hfs_freqs = [f * au.Unit(str(u)) for f, u in zip(hfs_data['hfs_freq'][idx_cols], hfs_data['unit'][idx_cols])]
                # append to list
                hfs_freq_list.append(hfs_freqs)

            # loop over lines
            for jj in range(n_lines):

                # line name from database
                line_name = lines_data["line_name"][jj]

                # create hyperfine structure mask for respective lines
                if line_name in lines_hfs:

                    print(f'{"[INFO]":<10}', f'Creating hyperfine structure mask for {line_name}.')

                    # get reference rest frequency and hyperfine frequencies
                    idx_line = lines_hfs.index(line_name)
                    hfs_reffreq = line_restfreq_list[idx_line]
                    hfs_freq = hfs_freq_list[idx_line]

                    # channel width
                    v_ch = this_data.meta["SPEC_DELTAV"]

                    # initionalise master finestructure mask (combination of all shifted masks)
                    mask_hfs = np.copy(mask)

                    # iterate over hyperfine frequencies
                    for freq, restfreq in zip(hfs_freq, hfs_reffreq):
                        
                        # compute velocity shift from frequency
                        freq_to_vel = au.doppler_radio(restfreq)
                        v_shift = freq.to(au.km/au.s, equivalencies=freq_to_vel) 
                        v_ch = v_ch.to(au.km/au.s)

                        # compute hyperfine velocity shift in amounts of channels
                        shift_ch = int(np.rint(v_shift.value / v_ch.value))

                        # shift mask to hyperfine frequency
                        mask_shift = np.empty_like(mask, dtype=float)
                        mask_shift[:] = 0
                        if shift_ch > 0:
                            mask_shift[:, shift_ch:] = mask[:, :-shift_ch]
                        elif shift_ch < 0:
                            mask_shift[:, :shift_ch] = mask[:, -shift_ch:]
                        else:
                            mask_shift = mask.copy()

                        # add to master mask
                        mask_hfs[mask_shift == 1] = 1

                    # assign units
                    mask_hfs *= au.dimensionless_unscaled

                    # save mask to database
                    this_data[f'SPEC_MASK_{line_name.upper()}'] = Column(mask_hfs, unit=au.dimensionless_unscaled, description=f'Velocity-integration mask for {line_name.upper()}')


        #store the mask in the PyStructure
        this_data["SPEC_MASK"] = Column(mask, unit=au.dimensionless_unscaled, description='Velocity-integration mask (used for integrated products)')
        #this_data["INT_VAL_VSHUFF"] = ref_line_vmean #JdB: remove, not needed in final product

        print(f'{"[INFO]":<10}', 'Done with mask. Computing moments.')

        #-------------------------------------------------------------------
        # Load noise mask (if noise velocity ranges were specified)
        #-------------------------------------------------------------------
        noise_mask = None
        if use_noise_vel_ranges and 'SPEC_NOISE_MASK' in this_data.keys():
            noise_mask = this_data['SPEC_NOISE_MASK']
            print(f'{"[INFO]":<10}', 'Using explicit noise velocity window(s) for RMS estimation.')

        #-------------------------------------------------------------------
        # Apply the mask to all lines and shuffle them
        #-------------------------------------------------------------------
        n_chan_new = 200 # LN: not used (remove?)
       
        for jj in range(n_lines):
            line_name = lines_data["line_name"][jj]

            if not 'SPEC_'+line_name.upper() in this_data.keys():
                print(f'{"[ERROR]":<10}', f'Tag for line {line_name.upper()} not found. Proceeding.')
                continue
            this_spec = this_data[f'SPEC_{line_name.upper()}']
            if np.nansum(this_spec, axis = None)==0:
                print(f'{"[ERROR]":<10}', f'Line {line_name.upper()} appears empty. Skipping.')
                continue

            dim_sz = np.shape(this_spec)
            n_pts = dim_sz[0]
            n_chan = dim_sz[1]
            this_v0 = this_data.meta["SPEC_VCHAN0"]
            this_deltav = this_data.meta["SPEC_DELTAV"]
            this_crpix = this_data.meta["SPEC_CRPIX"]
            
            this_vaxis = (this_v0 + (np.arange(n_chan)-(this_crpix-1))*this_deltav).to(au.km/au.s) #to km/s
            this_data["SPEC_VAXIS"] = Column(np.array([this_vaxis]*n_pts), unit=au.km/au.s, description='Velocity axis')

            # LN: Why do we have to do this? Should't the velocity axis of all lines be the same?
            # shuffled_mask = shuffle(spec = mask, \
            #                         vaxis = ref_line_vaxis,\
            #                         zero = 0.0,\
            #                         new_vaxis = this_vaxis, \
            #                         interp = 0)
                        
            # compute moment_maps
            # mom_maps = get_mom_maps(this_spec, shuffled_mask,this_vaxis, mom_calc)
            if use_hfs_lines:
                if line_name in lines_hfs:
                    mask_hfs = this_data[f'SPEC_MASK_{line_name.upper()}']* au.Unit(1)
                    mom_maps = get_mom_maps(this_spec, mask_hfs, this_vaxis, mom_calc, noise_mask=noise_mask)
                else:
                    mom_maps = get_mom_maps(this_spec, mask, this_vaxis, mom_calc, noise_mask=noise_mask)
            else:
                mom_maps = get_mom_maps(this_spec, mask, this_vaxis, mom_calc, noise_mask=noise_mask)

            # Save in structure
            line_desc = lines_data["line_desc"][jj]
            if lines_data["band_ext"].isnull()[jj]:
                
                tag_ii = "MOM0_"+line_name.upper()
                tag_uc = "EMOM0_" + line_name.upper()
                
                tag_tpeak = "TPEAK_" + line_name.upper()
                tag_rms = "RMS_" + line_name.upper()
                
                tag_mom1 = "MOM1_" + line_name.upper()
                tag_mom1_err = "EMOM1_" + line_name.upper()
                
                #Note that Mom2 corresponds to a FWHM
                tag_mom2 = "MOM2_" + line_name.upper()
                tag_mom2_err = "EMOM2_" + line_name.upper()
                
                tag_ew = "EW_" + line_name.upper()
                tag_ew_err = "EEW_" + line_name.upper()
                
                # store the different calculations
                this_data[tag_ii] = Column(mom_maps["mom0"], description=f'{line_desc} integrated intensity (moment-0)')
                this_data[tag_uc] = Column(mom_maps["mom0_err"], description=f'Propagated statistical error {line_desc} integrated intensity (moment-0)')
                this_data[tag_tpeak] =Column( mom_maps["tpeak"], description=f'{line_desc} peak brightness temperature')
                this_data[tag_rms] = Column(mom_maps["rms"], description=f'Statistical error {line_desc} peak brightness temperature')
                this_data[tag_mom1] = Column(mom_maps["mom1"], description=f'{line_desc} mean velocity (moment-1)')
                this_data[tag_mom1_err] = Column(mom_maps["mom1_err"], description=f'Propagated statistical error {line_desc} mean velocity (moment-1)')
                this_data[tag_mom2] = Column(mom_maps["mom2"], description=f'{line_desc} velocity dispersion (moment-2; {mom_calc[2]} definition)')
                this_data[tag_mom2_err] = Column(mom_maps["mom2_err"], description=f'Propagated statistical error {line_desc} velocity dispersion (moment-2; {mom_calc[2]} definition)')
                this_data[tag_ew] = Column(mom_maps["ew"], description=f'{line_desc} equivalent width (Gaussian approx)')
                this_data[tag_ew_err] = Column(mom_maps["ew_err"], description=f'Propagated statistical error {line_desc} equivalent width (Gaussian approx)')
                
                print(f'{"[INFO]":<10}', f'Done with line {lines_data["line_name"][jj]}.')
            
            else:
                print(f'{"[INFO]":<10}', f'Intensity Map for {lines_data["line_name"][jj]} already provided. Skipping.')

            #Shuffle the line
            #;- DC modify 02 march 2017: define a reference velocity axis
            #;-   this_deltav varies from dataset to dataset (fixing bug for inverted CO21 vaxis)
            cdelt = shuff_axis[1]* au.m / au.s
            naxis_shuff = int(shuff_axis[0])
            new_vaxis = cdelt * (np.arange(naxis_shuff)-naxis_shuff/2)
            new_vaxis=new_vaxis.to(au.km / au.s) #to km/s

            shuffled_line = shuffle(spec = this_spec,\
                                    vaxis = this_vaxis,\
                                    zero = ref_line_vmean,
                                    new_vaxis = new_vaxis,\
                                    interp = 0)

            tag_i = "SPEC_SHUFF" + line_name.upper()
            tag_v0 = "SPEC_VCHAN0_SHUFF"
            tag_deltav = "SPEC_DELTAV_SHUFF"


            this_data[tag_i] = Column(shuffled_line, unit=this_spec.unit,description=f'Shuffled {line_desc} brightness temperature')
            this_data.meta[tag_v0] = new_vaxis[0]
            this_data.meta[tag_deltav] = (new_vaxis[1] - new_vaxis[0])

            this_data["SPEC_VAXISSHUFF"] = Column(np.array([new_vaxis]*n_pts), unit=au.km/au.s, description="Shuffled Velocity Axis")
        

        this_data.meta["SPEC_CRPIX_SHUFF"] = 1
        this_data.write(fname[ii], format='ascii.ecsv', overwrite=True)

        print(f'{"[INFO]":<10}', f'Done with moments for {this_source}.')


        # /__
