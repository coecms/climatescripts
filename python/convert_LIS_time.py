#please use with care this has been put together quickly
#I'll be working on a better version to add to a github for anyone who needs to do the same
# Paola Petrelli paolap@utas.edu.au  20/10/2016
# This python script convertes the time variable of a LIS variable from the original format
# absolute time as double: YYYYMMDDHHmm with missing attributes
# to a CF compliant relative time : minutes since reference date 
# for details on how to run the script:
# python convert_LIS_time.py -h/--help
# This script works with both python2 and python3

from __future__ import print_function

import netCDF4 as nc4
import datetime as dt
import argparse
import sys
from glob import glob

def parse_input():
    ''' Parse input arguments '''
    parser = argparse.ArgumentParser(description=r'''Converts the time variable 
             in LIS output netcdf files from the original format 
             absolute time as double: YYYYMMDDHHmm with missing attributes
             to a CF compliant relative time : minutes since reference date 
            Takes as arguments the directory containing the input files 
            example to select two variables:
            python -d /short/w35/LIS_output
            It assumes the LIS files are named as 
            arguments are optional.
            The new variable time is written in the same input file, if you want
            to keep the original files unchanged, copy the files before running 
            the script.''',formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d','--input', type=str, nargs=1, help='input files directory', required=True)
    parser.add_argument('-s','--start_date', type=str, nargs=1, 
                        help='reference date to use in time units, format is yyyy-mm-dd',required=True)
    parser.add_argument('-fn','--filename', type=str, nargs=1, 
                        help='string matching filenames, this can be exactly one file or contain a wildcard "*" to' 
                        +'indicate more than one file, if not present "CABLE.*.nc" will be used instead',
                        required=False, default="CABLE.*.nc")
    return vars(parser.parse_args())


def time_to_date(value):
    ''' convert from time value to a datetime date '''
    val_str=str(int(value))
    yr=int(val_str[0:4])
    mn=int(val_str[4:6])
    dd=int(val_str[6:8])
    hh=int(val_str[8:10])
    mm=int(val_str[10:])
    return dt.datetime(yr,mn,dd,hh,mm)


def main():
    # assign input arguments
    args = parse_input()
    input_dir = args.pop("input")[0]
    start_date = args.pop("start_date")[0]
    filename = args.pop("filename")
    print(type(filename))
    if type(filename) is list: filename=filename[0]
    print(input_dir,start_date,filename)
    # create a list of input files
    print(input_dir+"/"+filename)
    print(glob(input_dir+"/"+filename))
    infiles=glob(input_dir+"/"+filename)
    # print files that will be chnaged and ask user if he/she wants to proceed
    print("Script will change time variable in the following files")
    for fname in infiles:
        print(fname)
    if sys.version_info < ( 3, 0 ):
        request=raw_input("Proceed? Y/N \n")
    else:
        request=input("Proceed? Y/N \n")
    if request == "Y":
        pass 
    else:
        sys.exit()
    for fname in infiles:
        # read the file and try to amend it
        nc = nc4.Dataset(fname,'r+')
        # read time variable from file
        time = nc.variables["time"]
        # define attributes for time variable
        time.standard_name = "time"
        time.units = "minutes since " + start_date + " 00:00"
        time.calendar = "standard"
        # read time values and convert them to new format 
        time_val = time[:]
        dates = []
        # first convert values to a "datetime" date format
        for val in time_val:
            dates.append(time_to_date(val))
        # then convert dates to numeric values
        time[:] = nc4.date2num(dates, units=time.units, calendar=time.calendar)
        print( time)
        nc.close()

if __name__ == "__main__":
    main()

        

