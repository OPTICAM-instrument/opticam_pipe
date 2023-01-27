import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]



def snr(rate,bkg,time,npix,rn,gain,dark=0.0,binning=1):
    source = rate * time
    shot_noise = rate * time
    sky_noise = bkg * npix * time*binning
    ro_noise = (rn**2 + (gain/2.0)**2 * npix*binning)
    dark_noise = dark * npix * time*binning
    
    return source / np.sqrt(shot_noise + sky_noise + ro_noise + dark_noise)

def snr_all(rate,bkg,time,npix,rn,gain,dark=0.0,binning=1):
    source = rate * time
    shot_noise = rate * time
    sky_noise = bkg * npix * time * binning
    ro_noise = (rn**2 + (gain/2.0)**2 * npix)*binning
    dark_noise = dark * npix * time * binning
    
    return source / np.sqrt(shot_noise ), source / np.sqrt(sky_noise),source / np.sqrt(ro_noise), source / np.sqrt(dark_noise)



import os
def rename_folder(folder):
    count = 0
    # count increase by 1 in each iteration
    # iterate all files from a directory
    for file_name in os.listdir(folder):
        # Construct old file name
        if not ('.fit' in file_name[-4:] or 'fits' in file_name[-4:]):
            continue
        print(file_name)
        name, extension = file_name.split('.')
        source = folder + file_name
        if 'u' in name or 'g' in name:
            ch = 'C1'
            if 'u' in name:
                key = 'u'
            else:
                key = 'g'
        if 'r' in name:
            ch = 'C2'
            key = 'r'
        if 'i' in name or 'z' in name:
            ch = 'C3'
            if 'i' in name:
                key = 'i'
            else:
                key = 'z'
        
        print(key)
        print(ch)
        print(name.split(key))
        # Adding the count to the new file name and extension
        new_name = name.split(key)[0]+ch+key+name.split(key)[1]+'.fits'
        destination = folder + new_name
    
        # Renaming the file
        os.rename(source, destination)
        count += 1
    print(f'All Files Renamed, \n Count:{count}')
    
    #print('New Names are')
    # verify the result
    #res = os.listdir(folder)
    #print(res)