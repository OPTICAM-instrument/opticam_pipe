import glob
import logging
import os
import sys
from datetime import datetime

import pandas as pd
from .misc import rename_folder
import sep


class Reduction:
    """
    Class for reduction and perform photometry for the
    data from the OPTICam instruments.

    Parameters
    ----------
    workdir : str
        Path to the working directory where the data products will be stored.

    rawdir : str
        Path to the directory where the raw data is stored.

    name : str
        Name of the target.

    camID: str
        This will be the camera id and will be used in the glob.glob rule.
        Default: 'C1'
        Default rule: '*C1*.fits'

    instrument: str
        This will be the instrument name and will be used to load the properties of the camera from the XX-properties.csv.


    ... more to come
    """

    def __init__(self, workdir, rawdir, name=None, camID='C1', instrument='MX'):

        # this should go here to do not create a conflict with the loops
        self.name = name

        # check if workdir ends in / and if not add it
        if workdir[-1] != '/':
            workdir += '/'
        # c check if workdir start with ~ and if yes replace it with the home directory
        if workdir[0] == '~':
            workdir = os.path.expanduser(workdir)

        self.workdir = workdir

        # check if workdir exists and if not create it after asking user
        if not os.path.exists(self.workdir):
            print('The directory {} does not exist.'.format(self.workdir))
            print('Do you want to create it? [y/n]')
            answer = input()
            if answer == 'y':
                # create the subdirectories needed
                os.makedirs(self.workdir)
                # create the logs directory
                os.makedirs(os.path.join(self.workdir, 'logs'))
                # log the creation of the directory
                print('Directory created.')

            else:
                print('Please enter a valid directory path.')
                sys.exit()

        # check if rawdir ends in / and if not add it
        if rawdir[-1] != '/':
            rawdir += '/'
        # check if rawdir start with ~ and if yes replace it with the home directory
        if rawdir[0] == '~':
            rawdir = os.path.expanduser(rawdir)
        # check if rawdir exists and if not exit
        if not os.path.exists(rawdir):
            print('The directory {} does not exist.'.format(rawdir))
            print('Please enter a valid directory path.')
            sys.exit()
        else:
            self.rawdir = rawdir

        # if the name of the target is none request it from the user
        if self.name is None:
            print('Please enter the name of the target.')
            self.name = input()
            print('Target name: {}'.format(self.name))
        else:
            print('Target name: {}'.format(self.name))

        # check if name ends in / if not add it
        if self.name[-1] != '/':
            self.name += '/'

        # create the subdirectories for target in the working directory
        self.targetdir = os.path.join(self.workdir, self.name)
        print('Target directory: {}'.format(self.targetdir))
        # check if the target directory exists and if not create it
        if not os.path.exists(self.targetdir):
            os.makedirs(self.targetdir)
            print('Target directory created.')
        else:
            print('Target directory already exists.')

        self.photdir = os.path.join(self.workdir, 'photometry')
        self.caldir = os.path.join(self.workdir, 'calibrated')
        self.logdir = os.path.join(self.workdir, 'logs')
        self.logfile = os.path.join(self.logdir, 'opticam.log')
        # check if the log file exist and if so request if the user wants to delete it, the default option is no
        if os.path.exists(self.logfile):
            print('The log file {} already exists.'.format(self.logfile))
            print('Do you want to delete it? [y/N]')
            answer = input()
            if answer == 'y':
                os.remove(self.logfile)
                print('Log file deleted.')
            else:
                pass
        self.log = logging.getLogger(__name__)
        self.log.setLevel(logging.DEBUG)
        self.log.addHandler(logging.FileHandler(self.logfile))
        self.log.addHandler(logging.StreamHandler())
        self.log.info('========== OPTICam pipeline ==========')
        # log current time and date
        self.log.info('Current date and time: {}'.format(datetime.now()))
        self.log.info('Starting reduction')
        self.log.info('Working directory: {}'.format(self.workdir))
        self.log.info('Target directory: {}'.format(self.targetdir))
        self.log.info('Raw data directory: {}'.format(self.rawdir))
        self.log.info('Calibrated data directory: {}'.format(self.caldir))
        self.log.info('Photometry data directory: {}'.format(self.photdir))
        self.log.info('Log file: {}'.format(self.logfile))
        self.log.info('Log level: {}'.format(logging.getLevelName(self.log.getEffectiveLevel())))
        self.log.info('Log level: {}'.format(logging.getLevelName(self.log.getEffectiveLevel())))

        # check if the names of the files are valid
        valid_names = ['C1', 'C2', 'C3']

        # check if the camera ID is valid
        if camID not in valid_names:
            print('Please enter a valid camera ID.')
            sys.exit()
        # report the camera ID
        self.camID = camID
        self.log.info('Camera ID: {}'.format(self.camID))
        # check if the instrument is valid
        if instrument not in ['MX', 'ARG']:
            print('Please enter a valid instrument.')
            sys.exit()
        # report the instrument
        self.instrument = instrument
        self.log.info('Instrument: {}'.format(self.instrument))

        # load the pixel scale from the config file
        properties = pd.read_csv('opticam/data/' + self.instrument + '-properties.csv', skiprows=1)
        self.ccd_pixscale = properties[self.camID + '_scale'][0]
        self.log.info('Pixel scale: {} arsec/pix'.format(self.ccd_pixscale))

        # check if the raw files contain the keys valid_names
        total_valid = 0
        for name in valid_names:
            total_valid += len(glob.glob(os.path.join(self.rawdir, '*' + name + '*.fit*')))

        if total_valid == 0:
            print('No files found with the valid names.'.format(name))
            # print the number of fits files in the directory
            print(
                'Number of fits files in the directory: {}'.format(len(glob.glob(os.path.join(self.rawdir, '*.fit*')))))

            # ask the user if he wants to rename the files in the folder to the valid names
            print('Do you want to rename the files in the directory to the valid names? [y/n]')
            answer = input()
            if answer == 'y':
                # rename the files in the directory using rename_folder
                rename_folder(self.rawdir)
                # log the renaming of the files
                self.log.info('All files in the directory renamed to valid names.')
            else:
                # exit the program
                print('Please rename the files in the directory to valid names.')
                # log the reason for exiting
                self.log.info('No files found with the valid names.')
                sys.exit()
        else:
            # log the number of files found
            self.log.info('Total number of valid found for all the cameras: {}'.format(total_valid))

        # use glob to find the files in the raw directory that contain the camera ID
        self.rawfiles = glob.glob(os.path.join(self.rawdir, '*' + self.camID + '*.fit*'))
        # report the directory and the number of files found
        self.log.info('Raw data directory: {}'.format(self.rawdir))
        self.log.info('Target directory: {}'.format(self.targetdir))
        self.log.info("Number of files found for {}: {}".format(self.camID, len(self.rawfiles)))

    def extract(self):
        """
        Routine that uses sep (python version of Source Extractor) to perform
        photometry and create a catalogue of stars for each file.
        """
        pass

class Background(sep.Background):
    """
    Class to create a background object using sep.Background.
    more info in https://sep.readthedocs.io/
    """
    pass