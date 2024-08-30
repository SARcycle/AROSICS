import argparse
import datetime

import utils


class Parameters:
    def __init__(self):
        self.args = self.parse_parameters()
        self.check_parameters()

    def parse_parameters(self):
        parse = argparse.ArgumentParser(
            description="Program for coregistering satellite imagery to a reference image over Switzerland using AROSICS")
        parse.add_argument('-s', '--start', type=str, help='Start date in format YYYYMMDD', default=None)
        parse.add_argument('-e', '--end', type=str, help='End date in format YYYYMMDD', default=None)
        parse.add_argument('-u', '--upload', action='store_true',
                           help='Toggle for upload of dx/dy into Google Earth Engine')
        parse.add_argument('-m', '--mosaicing', action='store_true',
                           help='Toggle for mosaicing of tiles before calculation of dx/dy')
        parse.add_argument('-c', '--coregistration', action='store_true', help='Toggle for calculating dx/dy')

        args = parse.parse_args()
        self.__setattr__('start_date', args.start)
        self.__setattr__('end_date', args.end)
        self.__setattr__('upload', args.upload)
        self.__setattr__('mosaicing', args.mosaicing)
        self.__setattr__('coreg', args.coregistration)

    def check_parameters(self):
        if self.start_date is not None:
            try:
                self.__setattr__('start_date', utils.parse_date(self.start_date))
            except Exception as e:
                assert False, f'{e}'
        else:
            self.__setattr__('start_date',
                             (datetime.datetime.now() - datetime.timedelta(days=5)).replace(hour=0, minute=0, second=0,
                                                                                            microsecond=0))

        if self.end_date is not None:
            try:
                self.__setattr__('end_date',
                                 utils.parse_date(self.end_date) + datetime.timedelta(days=1) - datetime.timedelta(
                                     seconds=1))
            except Exception as e:
                assert False, f'{e}'
        else:
            self.__setattr__('end_date', datetime.datetime.now().replace(hour=23, minute=59, second=59, microsecond=0))

        if self.upload is None:
            self.__setattr__('upload', False)

        if self.mosaicing is None:
            self.__setattr__('mosaicing', False)

        if self.coreg is None:
            self.__setattr__('coreg', False)
