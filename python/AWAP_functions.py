# functions to querry AWAP data which lives at "/srv/ccrc/data02/z3236814/data/AWAP/DAILY/netcdf"
# JK (08/03/2013) - First attempt
# JFE (12/03/2013) - Made code "lighter", better calc of Julian day, fix indexing error
# JK(13/03/2013) - Clean functions and check
# plz email any bugs to Jatin.Kala.JK@gmail.com 

# import libs here
import sys # this is apparently needed to call exit()
import numpy as np # for array calcs
from datetime import date # use for dates handling later on ...
from calendar import monthrange # used to find number of days in any particular month, of any year
from scipy.io import netcdf_file # to open netcdf files etc

#--------------------------------------------------------------------------------------------------------
def get_data_any_day_awapV3(year,month,day,var='pre'): #set var with a default value
   """ This function get AWAP precip, tmax, or tmin for any day, of any given month, of any given year
       for which data exists. The data lives at "/srv/ccrc/data02/z3236814/data/AWAP/DAILY/netcdf".
       See Jason's CCRC Data wiki entry for more details. I do not give warnings if your date falls 
       outside data availability! plz email any bugs to Jatin.Kala.JK@gmail.com
       Function arguments:
       year - just the year, integer
       month - 1 -12, integer
       day - day of month, NOT julian day, integer
       var - must either be set to "pre", "tmax", or "tmin", string
    """ 
# lets make sure everything is in the right type just in case
   year = int(year)
   month = int(month)
   day  = int(day)
   var = str(var)
# JFE hereafter I've made the code lighter
# define base path for data
   fld='/srv/ccrc/data02/z3236814/data/AWAP/DAILY/netcdf/'
# define directories for each variable in a dictionary, var_path is added to fld later on to give full path
   var_path={'pre': 'Daily_calib_rainfall_analysis_V3/',\
            'tmin': 'Daily_minimum_temperature_analysis/',\
            'tmax': 'Daily_maximum_temperature_analysis/'}
# if var not defined properly, exit program   
   if var not in var_path:
      print "var must be either pre, tmax, or tmin, exiting prgram now"
      sys.exit()
   else:
# JFE this is a nice way to concatenate strings in Python
      data_path='%s%s' % (fld,var_path[var])
   
# now we have path, we open the right file. These are in "yearly" files, and idexed by the julian day, so we need to figure out the julian day first!
# JFE the date object is better
      Julian_day = date(year, month, day).timetuple()[-2] 

#   print(Julian_day)

# now open the right file, get complete name first
# JFE this is a nice way to concatenate strings in Python no. 2, with integers  
   file_name = '%s%s.%4i.nc' % (data_path,var,year)
# print(file_name)
   f_open = netcdf_file(file_name)
# get the data out
# JFE first index in Python is 0, so need to remove 1day, could do it before, or create a variable for the index
   val_return = f_open.variables[var].data[Julian_day-1,:,:]
# JFE please close to free memory
   f_open.close()
# return the data
   return(val_return)
#--------------------------------------------------------------------------------------------------------
def get_data_any_month_awapV3(year,month,var='pre'): #set var with a default value
   """ This function get AWAP precip, tmax, or tmin of any given month, of any given year
       for which data exists. For tmax or tmin, a monthly mean is calculated, for precip, a sum is done.
       The data lives at "/srv/ccrc/data02/z3236814/data/AWAP/DAILY/netcdf".
       See Jason's CCRC Data wiki entry for more details. I do not give warnings if your date falls 
       outside data availability! plz email any bugs to Jatin.Kala.JK@gmail.com
       Function arguments:
       year - just the year, integer
       month - 1 -12, integer
       var - must either be set to "pre", "tmax", or "tmin", string
    """
   year = int(year)
   month = int(month)
   var = str(var)

   fld='/srv/ccrc/data02/z3236814/data/AWAP/DAILY/netcdf/'
   var_path={'pre': 'Daily_calib_rainfall_analysis_V3/',\
            'tmin': 'Daily_minimum_temperature_analysis/',\
            'tmax': 'Daily_maximum_temperature_analysis/'}
   
   if var not in var_path:
      print "var must be either pre, tmax, or tmin, exiting prgram now"
      sys.exit()
   else:
      data_path='%s%s' % (fld,var_path[var])

# get julian day of 1st and last day of month
   Julian_day_start = date(year, month, 1).timetuple()[-2]
   number_days_in_month = monthrange(year,month)[1]
   Julian_day_end = date(year, month, number_days_in_month).timetuple()[-2]
#   print(Julian_day_start)
#   print(Julian_day_end)

# open file
   file_name = '%s%s.%4i.nc' % (data_path,var,year)
   f_open = netcdf_file(file_name)

# if tmax or tmin, take the mean, if pre, take sum
   if var == "tmax" or var == "tmin":
      val_return = np.mean(f_open.variables[var].data[Julian_day_start-1:Julian_day_end-1,:,:],0)
   else:
      val_retrun = np.sum(f_open.variables[var].data[Julian_day_start-1:Julian_day_end-1,:,:],0)

   f_open.close()
# return the value
#   print(val_return.shape)
   return(val_return)
#------------------------------------------------------------------------------------------------------
def get_data_any_season_awapV3(year,season,var='pre'):
   """ This function get AWAP precip, tmax, or tmin for any season of any given year
       for which data exists. The data lives at "/srv/ccrc/data02/z3236814/data/AWAP/DAILY/netcdf".
       See Jason's CCRC Data wiki entry for more details. I do not give warnings if your date falls 
       outside data availability! plz email any bugs to Jatin.Kala.JK@gmail.com
       This function simply usues get_data_any_month_awapV3 !! nothing special here
       Function arguments:
       year - just the year, integer
       season - either DJF, MMA, JJA, or SON, string
       var - must either be set to "pre", "tmax", or "tmin", string
    """
   # lets make sure everything is in the right type just in case
   year = int(year)
   season = str(season)
   var = str(var)
  
   # define a dictionary, which gives months and years for each season
   seas_months_years = {'DJF': np.array([[12,1,2],[year-1,year,year]]),\
           'MAM': np.array([[3,4,5],[year,year,year]]),\
           'JJA': np.array([[6,7,8],[year,year,year]]),\
           'SON': np.array([[9,10,11],[year,year,year]])}

   if season not in seas_months_years:
      print "season must be either DJF, MAM, JJA, or SON, exiting prgram now"
      sys.exit()
   else:
      months_years = seas_months_years[season]
      val_return = np.mean(np.array([get_data_any_month_awapV3(months_years[1,0],months_years[0,0],var), \
                                     get_data_any_month_awapV3(months_years[1,1],months_years[0,1],var), \
                                     get_data_any_month_awapV3(months_years[1,2],months_years[0,2],var)]),0) 
   print(val_return.shape)
   return(val_return) 
#-----------------------------------------------------------------------------------------------------
