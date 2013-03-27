#!/usr/bin/python
### Python script that uses all nc_var_tools.py
# L. Fita, CCRC - ARC CoE CSS, UNSW, Sydney, Australia
####### ###### ##### #### ### ## #
## g.e. # nc_var.py -f R1_CCRC_NARCliM_MOM_1950-2009_pracc.nc -o addvattrk -S 'calendar|standard|S' -v time 
## g.e. # nc_var.py -f R1_CCRC_NARCliM_MOM_1950-2009_pracc.nc -o out -S 1:-1 -v time 
## g.e. # nc_var.py -f R1_CCRC_NARCliM_MOM_1950-2009_pracc.nc -o mname -S rain -v pracc
## g.e. # nc_var.py -f R1_CCRC_NARCliM_MOM_1950-2009_pracc.nc -o addvattr -S 'comment|Lluis!Fita-123456' -v pracc
## g.e. # nc_var.py -f R1_CCRC_NARCliM_MOM_1950-2009_pracc.nc -o rmvattr -S 'comment' -v pracc
## g.e. # nc_var.py -f R1_CCRC_NARCliM_MOM_1950-2009_pracc.nc -o subsetyrs -v pracc -S 1979-2005
## g.e. # nc_var.py -f R1_CCRC_NARCliM_MOM_1950-2009_pracc.nc -o infvattrs -v pracc

from optparse import OptionParser
import numpy as np
from netCDF4 import Dataset as NetCDFFile
import os
import re
import nc_var_tools as ncvar

errormsg='nc_var.py: ERROR -- error -- ERROR -- error'
warnmsg='nc_var.py: WARNING -- warning -- WARNING -- warning'

### Options
##string_operation="operation to make: " + '\n' + " out, output values -S inidim1,[inidim2,...]:enddim1,[enddim2,...]"
helpvar="""operation to make: 
  addfattr, add a attributes to a variable from another file and -S [file]:[var]
  addfdim, add a new dimension from another file and -S [file]:[dimension]
  addfgattr, add global attribute from another fiel and -S [file]
  addfvar, add a new variable from another file and -S [file]:[variable]
  addgattr, add a new global attribute: addatr -S [attrname]|[attrvalue]
  addgattrk, add a new global attribute: addatr -S [attrname]|[attrvalue]|[kind(S (!, white spaces),I,R,D)]
  addref, add a new variable with dimension and attributes from an already existing 'variable ref' in the file and -S [variable ref]:[attr name]@[value]:[[attr2]@[value2], ...]:[value/file with values]
  mname, modify name -S newname
  addvattr, add a new attribute to any given variable: addvattr -S [attrname]|[attrvalue]
  addvattrk, add a new attribute to any given variable: addvattrk -S [attrname]|[attrvalue]|[kind(S (!, white spaces),I,R,D)]
  infgattrs, give the values of all the global attributes of the file
  infsinggattrs, give the value of a single global attribute of the file
  infsingvattrs, give the value of a single attribute of the variable
  infvars, give the names of all the variables of the file
  infvattrs, give the values of all the attributes of a variable
  mdname, modify dimension name -S olddname:newdname
  means, computes the spatial mean of the variable
  meant, computes the temporal mean of the variable -S [power](power of the polynomial fit)
  mname, modify name -S newname
  out, output values -S inidim1,[inidim2,...]:enddim1,[enddim2,...]
  rgattr, read a global attribute: rgattr -S [attrname]
  rvattr, read a variable attribute: rvattr -S [attrname]
  rmgattr, remove a global attribute: rmgattr -S [attrname]
  rmvariable, remove a variable: rmvariable
  rmvattr, remove an attribute to any given variable: rmvattr -S [attrname]
  subsetmns, retrieve a subset of values based on months: subsetmns -S [MM1]:[...:[MMn]] or [MMi]-[MMe] for a period (output as 'subsetmns.nc')
  subsetyrs, retrieve a subset of values based on years: subsetyrs -S [YYYY1]:[...:[YYYYn]] or [YYYYi]-[YYYYe] for a period (output as 'subsetyrs.nc')
  valmod, modifiy values of variable -S [modification]:
     sumc,[constant]: add [constant] to variables values
     subc,[constant]: substract [constant] to variables values
     mulc,[constant]: multipy by [constant] to variables values
     divc,[constant]: divide by [constant] to variables values
     lowthres,[threshold],[newvalue]: modify all values below [threshold] to [newvalue]
     upthres,[threshold],[newvalue]: modify all values above [threshold] to [newvalue]

'addfattr': fattradd(opts.varname, opts.values, opts.ncfile)
'addfdim': fdimadd(opts.values, opts.ncfile)
'addfgattr': fgaddattr(opts.values, opts.ncfile)
'addfvar': fvaradd(opts.values, opts.ncfile)
'addgattr': gaddattr(opts.values, opts.ncfile)
'addgattrk': gaddattrk(opts.values, opts.ncfile)
'addref': varaddref(opts.values, opts.ncfile, opts.varname)
'addvattr': varaddattr(opts.values, opts.ncfile, opts.varname)
'addvattrk': varaddattrk(opts.values, opts.ncfile, opts.varname)
'infgattrs': igattr(opts.ncfile)
'infsinggattrs': isgattr(opts.values, opts.ncfile)
'infsingvattrs': isgattr(opts.values, opts.ncfile, opts.varname)
'infvars': ivars(opts.ncfile)
'infvattrs': ivattr(opts.ncfile, opts.varname)
'mdname': chdimname(opts.values, opts.ncfile, opts.varname)
'means': spacemean(opts.ncfile, opts.varname)
'meant': timemean(opts.values, opts.ncfile, opts.varname)
'mname': chvarname(opts.values, opts.ncfile, opts.varname)
'rgattr': grattr(opts.values, opts.ncfile)
'rvattr': vrattr(opts.values, opts.ncfile, opts.varname)
'rmgattr': grmattr(opts.values, opts.ncfile)
'rmvariable': varrm(opts.ncfile, opts.varname)
'rmvattr': varrmattr(opts.values, opts.ncfile, opts.varname)
'subsetmns': submns(opts.values, opts.ncfile, opts.varname)
'subsetyrs': subyrs(opts.values, opts.ncfile, opts.varname)
'out': varout(opts.values, opts.ncfile, opts.varname)
'valmod': valmod(opts.values, opts.ncfile, opts.varname)
"""

parser = OptionParser()
parser.add_option("-f", "--netCDF_file", dest="ncfile", 
                  help="file to use", metavar="FILE")
parser.add_option("-o", "--operation", type='choice', dest="operation", choices=['addfattr', 'addfdim', 'addfgattr', 'addfvar', 'addgattr', 
    'addgattrk', 'addref', 'addvattr', 'addvattrk', 'infgattrs', 'infsinggattrs', 'infsingvattrs', 'infvars', 'infvattrs', 'mdname', 'means', 'meant', 'mname', 'out', 'rgattr', 'rvattr', 'rmvariable', 'rmvattr', 
   'rmgattr', 'subsetmns', 'subsetyrs', 'valmod'], 
                  help=helpvar, metavar="OPER")
parser.add_option("-S", "--valueS (when applicable)", dest="values", 
                  help="values to use according to the operation", metavar="VALUES")
parser.add_option("-v", "--variable", dest="varname",
                  help="variable to check", metavar="VAR")

(opts, args) = parser.parse_args()

#######    #######
## MAIN
    #######

####### ###### ##### #### ### ## #

varn=opts.varname
oper=opts.operation

if not os.path.isfile(opts.ncfile):
  print errormsg
  print '  File ' + opts.ncfile + ' does not exist !!'
  print errormsg
  quit()    

if oper == 'addfattr':
  ncvar.fattradd(opts.varname, opts.values, opts.ncfile)
elif oper == 'addfdim':
  ncvar.fdimadd(opts.values, opts.ncfile)
elif oper == 'addfgattr':
  ncvar.fgaddattr(opts.values, opts.ncfile)
elif oper == 'addfvar':
  ncvar.fvaradd(opts.values, opts.ncfile)
elif oper == 'addgattr':
  ncvar.gaddattr(opts.values, opts.ncfile)
elif oper == 'addgattrk':
  ncvar.gaddattrk(opts.values, opts.ncfile)
elif oper == 'addref':
  ncvar.varaddref(opts.values, opts.ncfile, opts.varname)
elif oper == 'addvattr':
  ncvar.varaddattr(opts.values, opts.ncfile, opts.varname)
elif oper == 'addvattrk':
  ncvar.varaddattrk(opts.values, opts.ncfile, opts.varname)
elif oper == 'infgattrs':
  ncvar.igattrs(opts.ncfile)
elif oper == 'infsinggattrs':
  ncvar.isgattrs(opts.values, opts.ncfile)
elif oper == 'infsingvattrs':
  ncvar.isvattrs(opts.values, opts.ncfile, opts.varname)
elif oper == 'infvars':
  ncvar.ivars(opts.ncfile)
elif oper == 'infvattrs':
  ncvar.ivattrs(opts.ncfile, opts.varname)
elif oper == 'mdname':
  ncvar.chdimname(opts.values, opts.ncfile, opts.varname)
elif oper == 'means':
  ncvar.spacemean(opts.ncfile, opts.varname)
elif oper == 'meant':
  ncvar.timemean(opts.values, opts.ncfile, opts.varname)
elif oper == 'mname':
  ncvar.chvarname(opts.values, opts.ncfile, opts.varname)
elif oper == 'rgattr':
  ncvar.grattr(opts.values, opts.ncfile)
elif oper == 'rvattr':
  ncvar.vrattr(opts.values, opts.ncfile, opts.varname)
elif oper == 'rmgattr':
  ncvar.grmattr(opts.values, opts.ncfile)
elif oper == 'rmvariable':
  ncvar.varrm(opts.ncfile, opts.varname)
elif oper == 'rmvattr':
  ncvar.varrmattr(opts.values, opts.ncfile, opts.varname)
elif oper == 'subsetmns':
  ncvar.submns(opts.values, opts.ncfile, opts.varname)
elif oper == 'subsetyrs':
  ncvar.subyrs(opts.values, opts.ncfile, opts.varname)
elif oper == 'out':
  ncvar.varout(opts.values, opts.ncfile, opts.varname)
elif oper == 'valmod':
  ncvar.valmod(opts.values, opts.ncfile, opts.varname)
else:
  print errormsg
  print '   The operation ' + oper + ' is not ready !!'
  print errormsg
  quit()
