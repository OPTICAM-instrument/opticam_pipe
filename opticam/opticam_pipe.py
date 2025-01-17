import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import glob
from pathlib import Path
import pandas as pd
import os
import aplpy
from astropy.time import Time
import astroalign as aa
import sys
from .misc import * #this is to sort the text using the numbers in it  
from astropy.nddata import CCDData
from astropy.table import Table

#%%%
class Reduction:
    '''
    Object that performs the aperture photometry
    with SExtractor and creates a list of stars.

    The `Reduction` class is used to compute the 
    differential photometry over a particular dataset

    Parameters
    ----------
    workdir : str, optional
        Working directory where the data and catalogues are stored

    catalogue : str, optional
        Directory of all the stars catalogue as measured by SExtractor

    name : str, optional
        Name of the target

    rule : str, optional
        File rule to be used when collecting all the fits files. 
        Default '*.fits'

    Attributes
    ----------

    aper_size : float
        Aperture size used in the 

    airmass : float
        The airmass. This parameter is related to the altitude of the target.

    raw_data : data frame
        Dataframe that contains all the photometry from all stars in the field.

    all_stars : int, arrray
        Unique identifier of each star in the field

    epochs : int, array
        Array of epochs/images

    num_epochs : int
        Total number of epochs/images

    comp_factor : float, array
        Transmission factors at each epoch
    '''
#%%
    def __init__(self,workdir=None,rawdata = None,catalogue=None,
                name=None,rule='*.fits',config_fl_name=None, measurement_id=None, sizes=None,vrb=True):
        
        self.vrb = vrb
        if workdir is None: 
            self.workdir = './'
        else:
            self.workdir = workdir
            if not os.path.isdir(self.workdir):
                from pathlib import Path
                os.makedirs(self.workdir, exist_ok=True)
        if catalogue is None: 
            self.catalogue = 'catalogues/'
        else:
            self.catalogue = catalogue
        if rawdata is None: 
            self.rawdata = 'raw_data/'
        else:
            self.rawdata = rawdata
        if name is None:
            self.name = 'astro'
        else:
            self.name = name
        if config_fl_name == None:
            self.config_fl_name = 'default.sex'
        else:
            self.config_fl_name = config_fl_name
        if measurement_id is None:
            self.measurement_id = 'ISOCOR'
        elif measurement_id != 'ISO' and measurement_id != 'ISOCOR' and measurement_id != 'AUTO' and measurement_id != 'BEST' and measurement_id != 'APER' and measurement_id != 'PETRO':
            print('The inputed parameter does not correspond to any existing SExtractor parameters. Setting to default.')
            self.measurement_id = 'ISOCOR'
        else:
            self.measurement_id = measurement_id
            
        self.sizes = np.array([3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33])
        if isinstance(sizes,type(None)):
            print('No apperture size imputed, setting to default (16 pixels)')
            self.aper_ind = np.argwhere(self.sizes==16)[0]
        else:
            #here we implement several apertures
            self.aper_ind = [np.argwhere(np.array(self.sizes)==x)[0][0] for x in sizes]
            
        #if measurement_id == 'APER' and size is not None:
        #   
        #   if size in self.sizes:    
        #       self.aper_ind = self.sizes.index(size)
        #   else: #THIS NEED TO VE REVIEWED
        #       self.aper_ind = -1
        #       self.config_fl_name = self.config_fl_name.split('.')[0]+'_edit.sex'
        #       self.edit_sex_param(self.config_fl_name, ['PHOT_APERTURES'], [size])
                
        ###self.path_ref_list = False #inicializing the reference stars list
        
        self.rule = rule
        self.marker = '_C'+rule.split('C')[1][0]
        self.flns = self.get_files(self.rule)
        self._ROOT = os.path.abspath(os.path.dirname(__file__))
        self.path_ref_list = self.workdir+self.name+'_files/'+self.name+self.marker+'_ref_stars.csv'
#%%    

        #setting the pixelscale in the header
        if 'C1' in self.rule:
            self.ccd_pixscale = 0.1397
        elif 'C2' in self.rule:
            self.ccd_pixscale = 0.1406
        elif 'C3' in self.rule:
            self.ccd_pixscale =0.1661
        else:
            self.ccd_pixscale = 0.14
            
        
    def get_optimal_aperture(self):
        
        print('estimated FWHM arcseconds')
        print('mean:',np.mean(self.fwhm_image))
        print('median:',np.median(self.fwhm_image))
        print('min:',np.min(self.fwhm_image))
        print('max:',np.max(self.fwhm_image))
        print('std:',np.std(self.fwhm_image))
        print('\n\n')
        print('estimated FWHM pixels')
        print('mean:',np.mean(self.fwhm_image_pix))
        print('median:',np.median(self.fwhm_image_pix))
        print('min:',np.min(self.fwhm_image_pix))
        print('max:',np.max(self.fwhm_image_pix))
        print('std:',np.std(self.fwhm_image_pix))
        print('\n\n')
        print('binning:',np.unique(self.binning))
        print('ccd pixel scale:',self.ccd_pixscale)
        print('\n\n')        
        print('available apertures in pixels')
        print(self.sizes)
        
    def set_apertures(self,sizes):
        """
        Input:  
            sizes: array like of int
                aperture diameters in pixels 
        """
        print('updating aperture diameters (in pixels) to be saved')
        self.aper_ind = [np.argwhere(np.array(self.sizes)==x)[0][0] for x in sizes]
        
    def read_sex_param(self,fl_name):
        text = open(fl_name, 'r')
        lines = [line for line in text.readlines() if line.strip()]
        text.close()
        variables = []
        values = []
        for i in range(len(lines)):
            if lines[i][0] != '#':
                split = lines[i].split('\t')
                variables.append(split[0])
                if split[1] == '':
                    values.append(split[2])
                else:
                    values.append(split[1])
        d = {'Variables': variables, 'Values': values}
        dictionary = pd.DataFrame(data = d, dtype ='str')
        return dictionary
#%%
    def edit_sex_param(self, param, values, overwrite=False):
        
        fl_name = self.workdir+self.config_fl_name
        if not Path(fl_name).exists() or overwrite:
            sext_def= self._ROOT+'/sextractor_defaults/*'
            os.system('cp '+sext_def+' '+self.workdir)
            if self.vrb:
                print('generating default sextractor files')
            
            
            
        default = self.read_sex_param(fl_name)
        
        #self.confg_file= default
        
        for i,par in enumerate(param):
            default['Values'][default.Variables == par] = values[i]
        np.savetxt(fl_name,default.values,fmt='%s', delimiter='\t')
        
        if self.vrb: print('default params edited')
            
        #for i in range(len(param)):
        #    ss = (default['Variables'] == param[i])
        #    if param[i] == 'PHOT_APERTURES':
        #        default['Values'][ss] += ',' + str(values[i])
        #    else:
        #        default['Values'][ss] = values[i]
        #    if all(ss == False):
        #        d = {'Variables': [param[i]], 'Values': [values[i]]}
        #        d = pd.DataFrame(data = d, dtype ='str')
        #        default = pd.concat([default,d], ignore_index=True)
        #np.savetxt(fl_name,default.values,fmt='%s', delimiter='\t')
        return    
    
#%%%
    def get_files(self,rule):

        print('Looking in: ',self.workdir+self.rawdata+rule)
        self.flns = np.sort(glob.glob(self.workdir+self.rawdata+rule))

        if len(self.flns) == 0: 
            print('WARNING! >> No fits files detected')
        else:
            print('Found {} fits files.'.format(len(self.flns)))

        return self.flns

#%%%
    def sextractor(self):
        """
        Routine that uses SExtractor to perform
        aperture photometry and create a catalogue of 
        stars for each file.
        """
        current_dir = os.getcwd()
        
        os.chdir(self.workdir)
        
        self.binning = []
        self.fwhm_image = []
        self.fwhm_image_pix = []
        
        print(os.getcwd())
        fl_name_conf = self.config_fl_name
        if not Path(fl_name_conf).exists():
            sext_def= self._ROOT+'/sextractor_defaults/*'
            os.system('cp '+sext_def+' '+self.workdir)
            if self.vrb:
                print('generating default sextractor files') 
        else:
            if self.vrb: print('using existing sextractor files')
                
        if not os.path.isdir(self.catalogue):
            os.system('mkdir -p '+self.catalogue)

        for i,fln in enumerate(np.sort(glob.glob(self.rawdata+self.rule))):
            if fln[-4:] == "fits":
                cat_fln = fln.split(".fits")[0]+"_cat.fits"
            elif fln[-3:] == "fit":
                cat_fln = fln.split(".fit")[0]+"_cat.fits"
            exists = os.path.isfile(self.catalogue+ \
                    (cat_fln).split("/")[-1])
            
            if not exists:
                os.system("rm temp_sextractor_file.fits")

                hdul = fits.open(fln)
                #hotfix for broken files
                try: hdu1 = fits.PrimaryHDU(data=hdul[0].data)
                except: continue
                hdu1.header = hdul[0].header
                new_hdul = fits.HDUList([hdu1])
                new_hdul.writeto('temp_sextractor_file.fits', overwrite=True)

                hdul.close()
                os.system("chmod 777 temp_sextractor_file.fits")
                try:
                    gain = fits.getval(fln,"GAIN",0)
                except:
                    gain = 1.0
                    
                try:
                    self.binning.append(fits.getval(fln,"BINNING",0))
                    bnn = float(self.binning[-1].split('x')[0])
                except:
                    bnn = 1.
                    print("WARNING: BINNING NOT FOUND")

                sex_out = "sextractor temp_sextractor_file.fits  -c "+self.config_fl_name+" -CATALOG_NAME "+ \
                          cat_fln+" -GAIN "+str(gain)
                os.system(sex_out)
                print(cat_fln)
                
                #saving the estimated fwhm of the image
                msk = np.argwhere(fits.getdata(cat_fln).FWHM_IMAGE >0 ).T[0]
                PSF_FWHM_pix = np.median(fits.getdata(cat_fln).FWHM_IMAGE[msk])
                PSF_FWHM = PSF_FWHM_pix*self.ccd_pixscale * bnn
                
                self.fwhm_image.append(PSF_FWHM)
                self.fwhm_image_pix.append(PSF_FWHM_pix)
                ###
                
                #moving the cat file 
                print("mv "+cat_fln+\
                            " "+self.catalogue+".")
                os.system("mv "+cat_fln+\
                             " "+self.catalogue+".")

                print("{:4.0f} / {:4.0f} -- {}".format(i+1,len(self.flns),fln))
                
                
                
            else:
                print("{:4.0f} / {:4.0f} -- It exists!".format(i+1,len(self.flns)))
                
        os.chdir(current_dir)

#%%
    def creat_ref_list(self,number=0):
        '''
        Create reference star list

        number :: Default. First image of the list
        '''
        if not os.path.isdir(self.workdir+self.name+'_files/'):
            os.system('mkdir '+self.workdir+self.name+'_files/')
        fln = self.flns[number].split('/')[-1]

#        print(fln)
        if fln[-3:] == 'its':
            fl1 = self.workdir+self.catalogue+fln.split(".fits")[0]+"_cat.fits"
        else:
            fl1 = self.workdir+self.catalogue+fln.split(".fit")[0]+"_cat.fits"

        fl2 = self.workdir+self.rawdata+fln
        #print(fl1)

        data = fits.getdata(fl1)

        #print(data.columns)
        self.path_to_ref_fits = fl2

        fig = plt.figure(figsize=(14,14))
        gc = aplpy.FITSFigure(fl2,hdu=0,figure=fig)
        gc.show_grayscale(pmin=40,pmax=99,stretch="log",invert=True)

        gc.show_circles(data['X_IMAGE'], data['Y_IMAGE'], radius=13,color='g',lw=3)

        for i in range(data['X_IMAGE'].size):
            plt.text(data['X_IMAGE'][i]+10, data['Y_IMAGE'][i]+10,data['NUMBER'][i],fontsize=15,color='blue')
        
        plt.show()
        gc.savefig(self.workdir+self.name+'_files/'+self.name+self.marker+'_fov.pdf')
        df = pd.DataFrame(data=np.array([data['NUMBER'],
                                data['X_IMAGE'],
                                data['Y_IMAGE']]).T,columns=["id", "x","y"])
        
        self.path_ref_list = self.workdir+self.name+'_files/'+self.name+self.marker+'_ref_stars.csv'
        df.to_csv(self.path_ref_list,  index_label=False,index=False)

        self.ref_stars = df

#%%
    def get_position(self,num):

        ss = self.ref_stars.id == num

        self.tar_x = self.ref_stars.x.values[ss][0]
        self.tar_y = self.ref_stars.y.values[ss][0]
        print(self.tar_x,self.tar_y)
        print(self.name+ \
            ' position is: X= {:4.0f}, Y={:4.0f}'.format(self.tar_x,self.tar_y))

#%%

    def movie(self,target_id=None,clean_tmp=True):
        """
        Create a movie with all the images and the target cross matched. 
        This is based in the photometry method.
        
        target_id: index of the target in the reference image
                    
        clean_tmp: remove all the individual frames. Default = True
        """
        import gc as mpl
        self.photo_file = self.name+self.marker+'_photo' #+'_'+self.measurement_id
        apass = pd.read_csv(self.workdir+self.name+'_files/'+self.name+self.marker+'_ref_stars.csv',
            comment="#")
        apass.set_index('id')
        
        #if there is no target id we set it to 1 as default
        if not target_id:
            target_id = 1
                

        vrb = self.vrb 
        save_output = True


        PIX_EDGE = 30
        
        coo_apass = SkyCoord(apass['x']/1000.*u.deg, apass['y']/1000*u.deg)

        num_flns = len(self.flns)

        mov_fl = Path(self.workdir+self.name+'_files/'+self.name+self.marker+'_ref_stars.gif')
        print(mov_fl)

        if mov_fl.exists():
            print('Video file already exists')
            check_flag = True
        else:
            dd = {'flname': [], 'id_apass': [],'Filter': [],'MJD': [],
                 'flux': [],'flux_err': [],'mag': [],'mag_err': []
                 }

            sta = pd.DataFrame(data=dd)
            df3 = {}
            id3 = 0
            check_flag = False
        
        if check_flag: 
            print("File already exist")
            return 
    
        print("OPTICAM - Movie curve generator")
        
        ccd_pixscale = self.ccd_pixscale
        

        
        for i,flname in enumerate(self.flns[:]):
            k=i
            flnt = flname.split('/')[-1]
            if flnt.split('.')[-1] == 'fits':
                flnt2 = flnt[:-5]+'_cat.fits'
            elif flnt.split('.')[-1] == 'fit':
                flnt2 = flnt[:-4]+'_cat.fits'
                
            cat_flname = self.workdir+self.catalogue+flnt2
            #print(flname,cat_flname)
            if check_flag :
                if vrb: print(flname+" exists")
            else:
                print("Processing {:5.0f} / {:5.0f} : {}".format(i+1,num_flns,
                        flname.split('/')[-1]))

                filt = fits.getval(flname,"FILTER",0)
                #obj = fits.getval(flname,"OBJECT",0)
                exptime = fits.getval(flname,"EXPOSURE",0)
                try: mjd_t = fits.getval(flname,"GPSTIME",0)[:-5]
                except: mjd_t = fits.getval(flname,"UT",0)
                mjd_t = mjd_t.replace(' ', 'T')
                #hotfix 
                try: mjd = Time(mjd_t, format='fits', scale='utc').mjd
                except: #hotfix for new latest software version 
                    mjd_t =  fits.getval(flname,"DATE-OBS",0)+'T'+fits.getval(flname,"UT",0)
                    mjd = Time(mjd_t, format='fits', scale='utc').mjd
                airmass = fits.getval(flname,"AIRMASS",0)
                naxis1 = fits.getval(flname,"NAXIS1",0)
                naxis2 = fits.getval(flname,"NAXIS2",0)
                try: 
                    xbin= fits.getval(flname,"CCDXBIN",0)
                    ybin= fits.getval(flname,"CCDYBIN",0)
                    if xbin==ybin:
                        pixscale = ccd_pixscale * xbin
                    else:
                        pixscale = ccd_pixscale
                        print("Warning: different binning per axis")
                except:
                    pixscale= ccd_pixscale
                    print("Warning: Binning not found in the header, FWHM not trustable")
                   
                msk = np.argwhere(fits.getdata(cat_flname).FWHM_IMAGE >0 ).T[0]
                PSF_FWHM = np.median(fits.getdata(cat_flname).FWHM_IMAGE[msk])
                try:
                    seeing = fits.getval(flname,"L1FWHM",0)
                except:
                    seeing = PSF_FWHM*pixscale
                if seeing == "UNKNOWN": seeing = PSF_FWHM*pixscale
                if vrb: print("Seeing = {:7.3f} arcsec".format(seeing))
                if vrb: print("PSF FWHM = {:7.3f} arcsec".format(PSF_FWHM*pixscale))
                data = fits.getdata(cat_flname)


                #### Align images #####
                d_x,d_y = 0.0,0.0

                if i==0:
                    c_ref = np.array([data['X_IMAGE'],data['Y_IMAGE']]).T
                else:
                    c_tar = np.array([data['X_IMAGE'],data['Y_IMAGE']]).T
                    try:
                        p, (pos_img, pos_img_rot) = aa.find_transform(c_ref, c_tar)
                        d_x,d_y = p.translation[0],p.translation[1]
                        print("Translation: (x, y) = ({:.2f}, {:.2f})".format(*p.translation))
                    except: 
                        print('WARNING! >> List of matching triangles exhausted before an acceptable transformation was found?!?!')
                    

                coo_image = SkyCoord((data['X_IMAGE']-d_x)/1000*u.deg, 
                                     (data['Y_IMAGE']-d_y)/1000*u.deg)

                if vrb: print("Filter: {}".format(filt))

                idx_apass, d2d_apass, d3d_apass = coo_image.match_to_catalog_sky(coo_apass) #

                # Make mask due to separation
                ss = (d2d_apass.deg*1000 < 2)
                
                
                #print(idx_apass, d2d_apass, d3d_apass)
                #pass 
                std_mag = apass['x'][idx_apass]
                ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                measurement_id_list = ['ISO', 'ISOCOR', 'APER', 'AUTO', 'BEST', 'PETRO']
                flux_ISO = []
                flux_ISO_err = []
                mag_ISO = []
                mag_ISO_err = []
                
                flux_ISO.append(data['FLUX_ISO']/exptime)
                flux_ISO_err.append(data['FLUXERR_ISO']/exptime)
                mag_ISO.append(data['MAG_ISO'] + 2.5 * np.log10(exptime))
                mag_ISO_err.append(data['MAGERR_ISO'])
                
                flux_APER = []
                flux_APER_err = []
                mag_APER = []
                mag_APER_err = []
                
                flux_APER.append(data['FLUX_APER']/exptime)
                flux_APER_err.append(data['FLUXERR_APER']/exptime)
                mag_APER.append(data['MAG_APER'] + 2.5 * np.log10(exptime))
                mag_APER_err.append(data['MAGERR_APER'])

                pp = np.isfinite(np.array(mag_ISO)[0][ss]) & \
                     (data['X_IMAGE'][ss] > PIX_EDGE ) & \
                     (data['X_IMAGE'][ss] < naxis1 -PIX_EDGE)  & \
                     (data['Y_IMAGE'][ss] > PIX_EDGE ) & \
                     (data['Y_IMAGE'][ss] < naxis2 -PIX_EDGE )
                

                #creating individual plots for each image
                fig = plt.figure(figsize=(14,14))
                
                gc = aplpy.FITSFigure(flname,hdu=0,figure=fig,animated=True)
                gc.show_grayscale(pmin=40,pmax=99,stretch="log",invert=True)
                
                
                x_im, y_im = data['X_IMAGE'][ss][pp]+d_x, data['Y_IMAGE'][ss][pp]+d_y

                gc.show_circles(x_im, y_im , radius=13,color='g',lw=3,animated=True)
                
                #here we write the src index in each image corresponding to the reference image
                #we also identify the src index in the current image to substract the flux later. This is src_idx
                src_idx = None #reset this for each loop 
                for l in range(x_im.size):
                    plt.text(x_im[l]+6, y_im[l]+6,int(apass.id.iloc[std_mag[ss][pp].index.values[l]]),fontsize='large',animated=True)
                    
                    #rule 
                    if std_mag[ss][pp].index.values[l] == target_id-1:
                        gc.show_circles(x_im[l], y_im[l] , radius=14,color='r',lw=5,animated=True)
                        src_idx = l
                        
                    
                plt.text(10,10,flname.split('/')[-1],fontsize='large')    
                
                    
                
                
                #print()
                #print(self.aper_ind)
                #print(np.shape(mag_APER_err))
                #print(np.shape(pp),np.shape(pp))
                #print(mag_APER_err)
                #m_pl = mag_APER[0][:,self.aper_ind][ss][pp]
                #em_pl = mag_APER_err[0][:,self.aper_ind][ss][pp]
                m_pl = flux_APER[0][:,self.aper_ind][ss][pp]
                em_pl = flux_APER_err[0][:,self.aper_ind][ss][pp]
                #this is the size of the aperture in arcsec
                aper_size = pixscale * self.sizes[self.aper_ind]
                #plt.title('Airmass: {:.2f} SEEING: {:.2f} MAG aper: {:.2f} +/- {:.2f}, aper: {:.2f} arcsec'.format(airmass,seeing,m_pl[src_idx],em_pl[src_idx],aper_size))
                try:plt.title('Airmass: {:.2f} SEEING: {:.2f} Flux aper: {:.2e} +/- {:.2e}, aper: {:.2f} arcsec'.format(airmass,seeing,m_pl[src_idx],em_pl[src_idx],aper_size))
                except:
                    plt.title('ERROR with this frame')
                #return src_idx
                gc.savefig(self.workdir+self.name+'_files/fov{:1.0f}.jpg'.format(k))
                print(k)
                plt.close('all')
                plt.clf()
                mpl.collect() #this is to release the memory
                
                
                
        #creating the gif from jpgs   
        import glob
        from PIL import Image
        lis = [li.split('/')[-1] for li in glob.glob(self.workdir+self.name+'_files/*.jpg')]
        lis.sort(key=natural_keys)
        frames = [Image.open(self.workdir+self.name+'_files/'+image) for image in lis]
        frame_one = frames[0]
                                                          
        frame_one.save(self.workdir+self.name+'_files/'+self.name+self.marker+'_fov.gif', format="GIF", append_images=frames, save_all=True, fps=5, loop=0)
        
        #removing the temporal files
        if clean_tmp: 
            os.system('rm '+self.workdir+self.name+'_files/*.jpg')
        
                
                

        
    
    def photometry(self,PIX_EDGE = 30, vrb = None , save_output = True,save_standards = True,save_target = True):
        """
        Creates a single output file from all the catalogues. 
        Cross-matches the positions of each catalogue and assigns
        every star its unique identifier.
        
        PIX_EDGE: int, optional
            This avoid all the detections close to the edge of the CCD 
            default: ~4arsec ~30 pix 
        """
        self.photo_file = self.name+self.marker+'_photo' #+'_'+self.measurement_id
        apass = pd.read_csv(self.workdir+self.name+'_files/'+self.name+self.marker+'_ref_stars.csv',
            comment="#")
        apass.set_index('id')

        if vrb == None: vrb = self.vrb

        
        
        coo_apass = SkyCoord(apass['x']/1000.*u.deg, apass['y']/1000*u.deg)

        num_flns = len(self.flns)

        lco_log = Path(self.workdir+self.name+'_files/'+self.photo_file+'.csv')
        #print(lco_log)
        
        df3 = {}

        if lco_log.exists():
            print('Photometry file already exists')
            check_flag = True
        else:
            dd = {'flname': [], 'id_apass': [],'Filter': [],'MJD': [],
                 'flux': [],'flux_err': [],'mag': [],'mag_err': []
                 }

            sta = pd.DataFrame(data=dd)
            
            id3 = 0
            check_flag = False
        print("OPTICAM - Light curve generator")
        
        if 'C1' in self.rule:
            ccd_pixscale = 0.1397
        elif 'C2' in self.rule:
            ccd_pixscale = 0.1406
        elif 'C3' in self.rule:
            ccd_pixscale =0.1661
        else:
            ccd_pixscale = 0.14
        
        
        for i,flname in enumerate(self.flns[:]):
            flnt = flname.split('/')[-1]
            if flnt.split('.')[-1] == 'fits':
                flnt2 = flnt[:-5]+'_cat.fits'
            elif flnt.split('.')[-1] == 'fit':
                flnt2 = flnt[:-4]+'_cat.fits'
                
            cat_flname = self.workdir+self.catalogue+flnt2
            #print(flname,cat_flname)
            if check_flag :
                if vrb: print(flname+" exists")
            else:
                print("Processing {:5.0f} / {:5.0f} : {}".format(i+1,num_flns,
                        flname.split('/')[-1]))

                filt = fits.getval(flname,"FILTER",0)
                #obj = fits.getval(flname,"OBJECT",0)
                exptime = fits.getval(flname,"EXPOSURE",0)
                try: mjd_t = fits.getval(flname,"GPSTIME",0)[:-5]
                except: mjd_t = fits.getval(flname,"UT",0)
                mjd_t = mjd_t.replace(' ', 'T')
                #
                try: mjd = Time(mjd_t, format='fits', scale='utc').mjd
                except: #hotfix for new latest software version 
                    mjd_t =  fits.getval(flname,"DATE-OBS",0)+'T'+fits.getval(flname,"UT",0)
                    mjd = Time(mjd_t, format='fits', scale='utc').mjd
                airmass = fits.getval(flname,"AIRMASS",0)
                naxis1 = fits.getval(flname,"NAXIS1",0)
                naxis2 = fits.getval(flname,"NAXIS2",0)
                try: 
                    xbin= fits.getval(flname,"CCDXBIN",0)
                    ybin= fits.getval(flname,"CCDYBIN",0)
                    if xbin==ybin:
                        pixscale = ccd_pixscale * xbin
                    else:
                        pixscale = ccd_pixscale
                        print("Warning: different binning per axis")
                except:
                    pixscale= ccd_pixscale
                    print("Warning: Binning not found in the header, FWHM not trustable")
                    
                
                #creating a mask to elimitate 0 FWHM data
                #hotfix for bad data
                try: msk = np.argwhere(fits.getdata(cat_flname).FWHM_IMAGE >0 ).T[0]
                except: continue
                PSF_FWHM = np.median(fits.getdata(cat_flname).FWHM_IMAGE[msk])
                try:
                    seeing = fits.getval(flname,"L1FWHM",0)
                except:
                    seeing = PSF_FWHM*pixscale
                if seeing == "UNKNOWN": seeing = PSF_FWHM*pixscale
                if vrb: print("Seeing = {:7.3f} arcsec".format(seeing))
                if vrb: print("PSF FWHM = {:7.3f} arcsec".format(PSF_FWHM*pixscale))
                data = fits.getdata(cat_flname)
              


                #### Align images #####
                d_x,d_y = 0.0,0.0

                if i==0:
                    c_ref = np.array([data['X_IMAGE'],data['Y_IMAGE']]).T
                else:
                    c_tar = np.array([data['X_IMAGE'],data['Y_IMAGE']]).T
                    try:
                        p, (pos_img, pos_img_rot) = aa.find_transform(c_ref, c_tar)
                        d_x,d_y = p.translation[0],p.translation[1]
                        print("Translation: (x, y) = ({:.2f}, {:.2f})".format(*p.translation))
                    except: 
                        print('WARNING! >> List of matching triangles exhausted before an acceptable transformation was found?!?!')
                    

                coo_image = SkyCoord((data['X_IMAGE']-d_x)/1000*u.deg, 
                                     (data['Y_IMAGE']-d_y)/1000*u.deg)

                if vrb: print("Filter: {}".format(filt))

                idx_apass, d2d_apass, d3d_apass = coo_image.match_to_catalog_sky(coo_apass) #

                # Make mask due to separation
                ss = (d2d_apass.deg*1000 < 2)

                std_mag = apass['x'][idx_apass]
                ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                measurement_id_list = ['ISO', 'ISOCOR', 'APER', 'AUTO', 'BEST', 'PETRO']
                flux_ISO = []
                flux_ISO_err = []
                mag_ISO = []
                mag_ISO_err = []
                
                flux_ISO.append(data['FLUX_ISO'])
                flux_ISO_err.append(data['FLUXERR_ISO'])
                mag_ISO.append(data['MAG_ISO'] + 2.5 * np.log10(exptime))
                mag_ISO_err.append(data['MAGERR_ISO'])
                
                flux_ISOCOR = []
                flux_ISOCOR_err = []
                mag_ISOCOR = []
                mag_ISOCOR_err = []
                
                flux_ISOCOR.append(data['FLUX_ISOCOR'])
                flux_ISOCOR_err.append(data['FLUXERR_ISOCOR'])
                mag_ISOCOR.append(data['MAG_ISOCOR'] + 2.5 * np.log10(exptime))
                mag_ISOCOR_err.append(data['MAGERR_ISOCOR'])
                
                flux_APER = []
                flux_APER_err = []
                mag_APER = []
                mag_APER_err = []
                
                flux_APER.append(data['FLUX_APER'])
                flux_APER_err.append(data['FLUXERR_APER'])
                mag_APER.append(data['MAG_APER'] + 2.5 * np.log10(exptime))
                mag_APER_err.append(data['MAGERR_APER'])

                flux_APER = np.array(flux_APER)
                flux_APER_err = np.array(flux_APER_err)
                mag_APER = np.array(mag_APER)
                mag_APER_err = np.array(mag_APER_err)
            
                flux_AUTO = []
                flux_AUTO_err = []
                mag_AUTO = []
                mag_AUTO_err = []
                
                flux_AUTO.append(data['FLUX_AUTO'])
                flux_AUTO_err.append(data['FLUXERR_AUTO'])
                mag_AUTO.append(data['MAG_AUTO'] + 2.5 * np.log10(exptime))
                mag_AUTO_err.append(data['MAGERR_AUTO'])
                
                flux_BEST = []
                flux_BEST_err = []
                mag_BEST = []
                mag_BEST_err = []
                
                flux_BEST.append(data['FLUX_BEST'])
                flux_BEST_err.append(data['FLUXERR_BEST'])
                mag_BEST.append(data['MAG_BEST'] + 2.5 * np.log10(exptime))
                mag_BEST_err.append(data['MAGERR_BEST'])
                
                flux_PETRO = []
                flux_PETRO_err = []
                mag_PETRO = []
                mag_PETRO_err = []
                
                flux_PETRO.append(data['FLUX_PETRO'])
                flux_PETRO_err.append(data['FLUXERR_PETRO'])
                mag_PETRO.append(data['MAG_PETRO'] + 2.5 * np.log10(exptime))
                mag_PETRO_err.append(data['MAGERR_PETRO'])


                pp = np.isfinite(np.array(mag_ISO)[0][ss]) & \
                     (data['X_IMAGE'][ss] > PIX_EDGE ) & \
                     (data['X_IMAGE'][ss] < naxis1 -PIX_EDGE)  & \
                     (data['Y_IMAGE'][ss] > PIX_EDGE ) & \
                     (data['Y_IMAGE'][ss] < naxis2 -PIX_EDGE )
                

                if vrb: print("Number of Absolute detected stars {} \n ".format(pp.sum()))

                if ((pp.sum() >= 3)) & save_target:
                   if save_standards:
                        for jj in np.arange(pp.sum()):
                            df3[id3] = {'flname': flname, 'id_apass':apass.id.iloc[std_mag[ss][pp].index.values[jj]],
                                 'Filter': filt,'MJD': mjd+exptime/86400./2.,
                                 'epoch':i,
                                 'flux_ISO': flux_ISO[0][ss][pp][jj],
                                 'flux_err_ISO': flux_ISO_err[0][ss][pp][jj],
                                 'mag_ISO': mag_ISO[0][ss][pp][jj],
                                 'mag_err_ISO': mag_ISO_err[0][ss][pp][jj],
                                 'flux_ISOCOR': flux_ISOCOR[0][ss][pp][jj],
                                 'flux_err_ISOCOR': flux_ISOCOR_err[0][ss][pp][jj],
                                 'mag_ISOCOR': mag_ISOCOR[0][ss][pp][jj],
                                 'mag_err_ISOCOR': mag_ISOCOR_err[0][ss][pp][jj],
                                 'flux_AUTO': flux_AUTO[0][ss][pp][jj],
                                 'flux_err_AUTO': flux_AUTO_err[0][ss][pp][jj],
                                 'mag_AUTO': mag_AUTO[0][ss][pp][jj],
                                 'mag_err_AUTO': mag_AUTO_err[0][ss][pp][jj],
                                 'flux_BEST': flux_BEST[0][ss][pp][jj],
                                 'flux_err_BEST': flux_BEST_err[0][ss][pp][jj],
                                 'mag_BEST': mag_BEST[0][ss][pp][jj],
                                 'mag_err_BEST': mag_BEST_err[0][ss][pp][jj],
                                 'flux_PETRO': flux_PETRO[0][ss][pp][jj],
                                 'flux_err_PETRO': flux_PETRO_err[0][ss][pp][jj],
                                 'mag_PETRO': mag_PETRO[0][ss][pp][jj],
                                 'mag_err_PETRO': mag_PETRO_err[0][ss][pp][jj],
                                 'exptime': exptime,
                                 'airmass': airmass,
                                 'seeing':seeing
                                 }
                            
                            #here we save the different apertures:
                            for x, ap_ind in enumerate(self.aper_ind):
                                #print(x,ap_ind)
                                df3[id3][f'flux_APER_{x+1}']= flux_APER[0,:,ap_ind][ss][pp][jj]
                                df3[id3][f'flux_err_APER_{x+1}']= flux_APER_err[0,:,ap_ind][ss][pp][jj]
                                df3[id3][f'mag_APER_{x+1}']= mag_APER[0,:,ap_ind][ss][pp][jj]
                                df3[id3][f'mag_err_APER_{x+1}']= mag_APER_err[0,:,ap_ind][ss][pp][jj]
                                
                            #print()
                            id3 += 1
        #############################################################################################
        if (len(df3) >= 1) & save_target:
                sta = pd.DataFrame.from_dict(df3,"index")
                sta = sta.sort_values(by=['id_apass','epoch'])
                
                self.out_df = sta 
            
                    
                #here we update the targets that are saved but not detected in all the images
                if self.path_ref_list:
                    if vrb: 
                        print('adding number of detections to the id catalogue')
                        
                    
                    df_t = pd.read_csv(self.path_ref_list) #we load theref_star file
                    #n=np.zeros(len(df_t['id'])
                    #df_t['n']= np.zeros(len(df_t['id']))  #we wet all the observations to 0
                    for i in self.out_df.id_apass.unique(): #we sort the ids to match the reference file order
                        #n.append(len(self.out_df[(self.out_df.id_apass == i)]))
                               
                        n_i = len(self.out_df[(self.out_df.id_apass == i)]) #number of detections of the i istar in the catalogue
                               
                        #ID = np.argwhere(df_t['id'].values == i)
                        #df_t['detections'][ID]= n_i 
                        df_t.loc[df_t['id'].values == i, 'n'] = n_i
                    
                    
                                 
                    df_t.to_csv(self.path_ref_list)
                                 
                    if vrb: print('Done')
                    
                    #saving copying the headers of the reference image to the output files
                    ccd = CCDData.read(self.path_to_ref_fits, unit='count')
                    header_flag = True
                    sta.meta= ccd.meta[6:]
                    sta.meta['Camera'] = int(self.marker[-1])
                    
                    #saving number of pixels in the metadata 
                    for x, ap_ind in enumerate(self.aper_ind):
                        sta.meta[f'APER_{x+1}_d_pix'] = self.sizes[ap_ind]
                        #in arcsec
                        bnn = float(np.unique(self.binning)[-1].split('x')[0])
                        sta.meta[f'APER_{x+1}_d_pix'] = self.sizes[ap_ind] * bnn * self.ccd_pixscale
                        
                else: header_flag = False
                
                if save_output & save_standards:
                    sta.to_csv(self.workdir+self.name+'_files/'+self.photo_file+".csv")
                    sta.to_pickle(self.workdir+self.name+'_files/'+self.photo_file+".pkl")
                    
                    t = Table.from_pandas(sta)
                    t.meta = sta.meta
                    t.write(self.workdir+self.name+'_files/'+self.photo_file+".fits",overwrite=True)
                    
                    print('Files saved in '+self.workdir+self.name+'_files/'+self.photo_file)
                    
                  
                

