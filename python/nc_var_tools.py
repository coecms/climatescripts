# netCDF utilities
# L. Fita, CCRC - ARC CoE CSS, UNSW, Sydney, Australia
####### ###### ##### #### ### ## #
## valmod: Function to modify the value of a variable
## rangedim: Gives the instruction to retrieve values from a dimension of a variable
## varaddref: Function to add a variable in an existing file copying characteristics from an existing one
## varout: Function when we want to output variable values
## chdimname: Changing the name of the dimension
## chvarname: Changing the name of the variable
## searchInlist: Function to search a value within a list
## set_attribute: Sets a value of an attribute of a netCDF variable. Removes previous attribute value if exists
## set_attributek: Sets a value of an attribute of a netCDF variable with a kind. Removes previous attribute value if exists
## gaddattr: Add a global attribute to a netCDF. Removes previous attribute if it exist
## gaddattrk: Add a global attribute to a netCDF caring about the type. Removes previous attribute if it exist
## varaddattr: Add an attribute to a variable. Removes previous attribute if it exists
## varaddattrk: Add an attribute to a variable caring about the type
## varrmattr: Removing an attribute from a variable
## grmattr: Removing a global attribute
## fvaradd: Adding variable (and all its attributes and dimensions) from a reference file to a file
## fdimadd: Adding dimension from another reference file 
## fattradd: Adding attributes from a reference file
## fgaddattr: Adding global attributes from a reference file
## varrm: Removing a variable from a file
## ivars: Give all the variable names of a file
## igattrs: Give all the global attributes of a file
## isgattrs: Give a single global attribute of a file and its type
## ivattrs: Give all the attributes of a variable and its type
## isvattrs: Give a single attribute of a variable
## grattr: Function to read a global atribute
## vrattr: Function to remove an atribute from a variable
## datetimeStr_datetime: Function to transform a string date ([YYYY]-[MM]-[DD]_[HH]:[MI]:[SS] format ) to a date object
## dateStr_date: Function to transform a string date ([YYYY]-[MM]-[DD] format) to a date object
## timeStr_time: Function to transform a string date ([HH]:[MI]:[SS] format) to a time object
## time_information: Function to provide information about variable time
## CFtimes_datetime: Provide date/time array from a file with a series of netCDF CF-compilant time variable
## variable_inf: Class to provide the information from a given variable
## subyrs: Function to retrieve a series of years from a file
## submns: Function to retrieve a series of months from a file
## statsValWeigthed: Weigthed Statistics class
## stats2Val: two variables Statistics class
## statsVal: Statistics class
## spacemean: Function to retrieve a space mean series from a multidimensional variable of a file
## timemean: Function to retrieve a time mean series from a multidimensional variable of a file
##

import numpy as np
from netCDF4 import Dataset as NetCDFFile
import os
import re

errormsg='nc_var_tools.py: ERROR -- error -- ERROR -- error'
warnmsg='nc_var_tool.py: WARNING -- warning -- WARNING -- warning'

def valmod(values, ncfile, varn):
    """ Function to modify the value of a variable
    values = modins,modval1,[modval2,...]
      modins = instruction: 
        'sumc', add [modval1]
        'subc', substraction [modval1]
        'mulc', multiply by [modval1] 
        'divc', divide by [modval1] 
        'lowthres': modify all values below [modval1] to [modval2]
        'upthres': modify all values above [modval1] to [modval2]
    ncfile = netCDF file name
    varn = name of the variable
    """
    import numpy as np
    from netCDF4 import Dataset as NetCDFFile
#  #  errormsg='ERROR -- error -- ERROR -- error'
    vals = values.split(',')
    modins = vals[0]
    modval = float(vals[1])

    if not os.path.isfile(ncfile):
      print errormsg
      print '   valmod: File "' + ncfile + '" does not exist !!'
      print errormsg
      quit(-1)    

    ncf = NetCDFFile(ncfile,'a')

    if ncf.dimensions.has_key('plev'):
      # removing pressure level from variable name
      varn = re.sub("\d+", "", varn) 

    if not ncf.variables.has_key(varn):
      print errormsg
      print '   valmod: File does not have variable "' + varn + '" !!!!'
      print errormsg
      ncf.close()
      quit(-1)

    var = ncf.variables[varn]
    varshape = var.shape
    Ndims = len(varshape)
    
    varshapevals = list(varshape)
    varVal = var[:]

    if modins == 'sumc':
      varVal[:] = varVal[:] + modval
    elif modins == 'subc':
      varVal[:] = varVal[:] - modval
    elif modins == 'mulc':
      varVal[:] = varVal[:] * modval
    elif modins == 'divc':
      varVal[:] = varVal[:] / modval
    elif modins == 'lowthres':
      varVal2 = np.where(varVal[:] < float(vals[1]), float(vals[2]), varVal[:])
      varVal[:] = varVal2
    elif modins == 'upthres':
      varVal2 = np.where(varVal[:] > float(vals[1]), float(vals[2]), varVal[:])
      varVal[:] = varVal2
    else: 
      print errormsg
      print '  valmod: Operation to modify values ' + modins + ' is not defined !!!'
      print errormsg
      quit(-1)

    var[:] = varVal
    ncf.sync()
    ncf.close()

def rangedim(end, shape):
    """Gives the instruction to retrieve values from a dimension of a variable
    >>> print rangedim(-1, 15)
    15
    """
    if end == -1:
      return shape
    else:
      return end

def varaddref(values, ncfile, varn):
  """ Function to add a variable in an existing file copying characteristics from an existing one
  values = [variable ref]:[attr name]@[value][:[attr2]@[value2], ...]:[value/file with values] add a new variable [varn] 
    with dimension and attributes from an already existing [variable ref] with attributes [[attr name]@[value][:[attr2]@[value2], ...]] 
    in the file [netcdf] and value [value/file with valules]
  netcdf = netCDF file name
  varn = new variable name
  """

  varvalues = values.split(':')

  Nvarvalues = len(varvalues)
  varprev = varvalues[0]
  newattrs = {}
  for iattr in range(Nvarvalues - 2):
    attrv = varvalues[iattr+1]
    newattrs[attrv.split('@')[0]] = attrv.split('@')[1]

  ncf = NetCDFFile(ncfile,'a')

  if ncf.variables.has_key(varn):
    print errormsg
    print '   varaddref: File already has the varible ' + varn + ' !!!'
    print errormsg
    ncf.close()
    quit(-1)

  if not ncf.variables.has_key(varprev):
    print errormsg
    print '    varaddref: File does not have variable ' + varprev + ' !!!!'
    print errormsg
    ncf.close()
    quit(-1)

  varp = ncf.variables[varprev]
  varpVal = varp[:]
  varpshape = varp.shape
  Nvarpshape = len(varpshape)
  varpattrs = varp.ncattrs()
  varptype = varp.dtype
  varpdims = varp.dimensions
  varpndims = varp.ndim

##  print '  shape of the variable used as reference: ',varpshape
#
# variable characteristics 
  if len(varpshape) == 4:
    dimpt = varpshape[0]
    dimpz = varpshape[1]
    dimpy = varpshape[2]
    dimpx = varpshape[3]
  elif len(varpshape) == 3:
    dimpt = varpshape[0]
    dimpy = varpshape[1]
    dimpx = varpshape[2]
  elif len(varpshape) == 3:
    dimpy = varpshape[0]
    dimpx = varpshape[1]

  newvar = ncf.createVariable(varn, varptype, dimensions=varpdims)
  newvar = ncf.variables[varn]
#
# Values
  varplen=1
  for idim in range(varpndims):
    varplen=varplen*varpshape[idim]

  if varptype == 'float' or varptype == 'float32' or varptype == 'float64' or vartype == np.float(1.) or vartype == np.float32(1.):
    newvals = np.reshape(np.arange(varplen), varpshape)*1.
  else:
    newvals = np.reshape(np.arange(varplen), varpshape)

  if not os.path.isfile(varvalues[Nvarvalues-1]):
##    print '  Using constant value'
    newvals[:] = float(varvalues[Nvarvalues-1])
  else:
##    print '  Using 2-D values from ' + varvalues[Nvarvalues-1]
    asciif = open(varvalues[Nvarvalues-1], 'r')

    fdimx = len(asciif.readlines())
    asciif.seek(0)

    if len(varpshape) == 4:
      if not fdimx == dimpx:
        print errormsg
        print '   varaddref: Provided file has dimx=', fdimx, ' and variable has dimx=', dimpx, ' !!!'
        print errormsg
        ncf.close()
        quit(-1)

      iline = 0
      idx = 0
      for fline in asciif:
        line = fline.replace('\n','')
        yvals = []
        for iyval in line.split('     '):
          yvals.append(float(iyval))

        if iline == 0:
          fdimy = len(yvals)
          if not fdimy == dimpy:
            print errormsg
            print '    varaddref: Provided file has dimy=', fdimy, ' and variable has dimy= ', dimpy, ' !!!'
            print errormsg
            ncf.close()
            quit(-1)
        for it in range(dimpt):
          for iz in range(dimpz):
            newvals[it,iz,:,idx] = yvals
  
        idx = idx+1
        iline = iline+1

    elif len(varpshape) == 3:
      if not fdimx == dimpx:
        print errormsg
        print '    varaddref: Provided file has dimx=', fdimx, ' and variable has dimx=', dimpx, ' !!!'
        print errormsg
        ncf.close()
        quit(-1)

      iline = 0
      idx = 0
      for fline in asciif:
        line = fline.replace('\n','')
        yvals = []
        for iyval in line.split('     '):
          yvals.append(float(iyval))

        if iline == 0:
          fdimy = len(yvals)
          if not fdimy == dimpy:
            print errormsg
            print '    varaddref: Provided file has dimy=', fdimy, ' and variable has dimy= ',dimpy, ' !!!'
            print errormsg
            ncf.close()
            quit(-1)
        for it in range(dimpt):
          newvals[it,:,idx] = yvals

        idx = idx+1
        iline = iline+1
    elif len(varpshape) == 2:
      if not fdimx == dimpx:
        print errormsg
        print '    varaddref: Provided file has dimx=', fdimx, ' and variable has dimx=', dimpx, ' !!!'
        print errormsg
        ncf.close()
        quit(-1)

      iline = 0
      idx = 0
      for fline in asciif:
        line = fline.replace('\n','')
        yvals = []
        for iyval in line.split('     '):
          yvals.append(float(iyval))

        if iline == 0:
          fdimy = len(yvals)
          if not fdimy == dimpy:
            print errormsg
            print '    varaddref: Provided file has dimy=', fdimy, ' and variable has dimy= ',dimpy, ' !!!'
            print errormsg
            ncf.close()
            quit(-1)

        newvals[:,idx] = yvals
        idx = idx+1
        iline = iline+1

  asciif.close()
  newvar[:] = newvals

#
# Attributes
  for iattr in range(len(varpattrs)):
    attrn = varpattrs[iattr]
    attrval = varp.getncattr(attrn)

    if attrn in newattrs:
      newattr = newvar.setncattr(attrn, newattrs[attrn])
    else:
      newattr = newvar.setncattr(attrn, attrval)

  ncf.sync()
  ncf.close()

def varout(values, ncfile, varn):
  """ Function when we want to output variable values
  values = [optsIrange]:[optsErange]
    [optsIrange]: val1,val2,...,valN inital value for the 'N' dimensions of the variable
    [optsErange]: val1,val2,...,valN ending value for the 'N' dimensions of the variable
    -1 for all range of values
  ncfile = netCDF file name
  varn = variable name
  """
  import numpy as np
  from netCDF4 import Dataset as NetCDFFile

  optsIrange = values.split(':')[0]
  optsErange = values.split(':')[1]

  print optsIrange
  print optsErange

  inirange=optsIrange.split(',')
  endrange=optsErange.split(',')

  irange = [int(val) for val in inirange]
  erange = [int(val) for val in endrange]

  if not len(irange) == len(erange):
    print errormsg
    print '    varout: Different number of values in each range!'
    print '    varout: initial range: ' + irange
    print '    varout: ending range: ' + erange
    print errormsg
    quit(-1)
  else:
    ndims=len(irange)
 
  if not os.path.isfile(ncfile):
    print errormsg
    print '    varout: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'r')

  if ncf.dimensions.has_key('plev'):
    # removing pressure level from variable name
    varn = re.sub("\d+", "", varn) 

  if not ncf.variables.has_key(varn):
    print errormsg
    print '    varout: File "' + ncfile + '" does not have variable "' + varn + '" !!!!'
    print errormsg
    ncf.close()
    quit(-1)

  var = ncf.variables[varn]

  varshape = var.shape
  Nvarshape = len(varshape)
 
  if not Nvarshape == ndims:
    print errormsg
    print '    varout: Provided number of values of the range ' + ndims + ' is different of the shape of the variable ' + Nvarshape + ' !!!'
    print errormsg
    ncf.close()
    quit(-1)

  if Nvarshape == 1:
    varValrange = var[irange[0]:rangedim(erange[0], varshape[0])]
  elif Nvarshape == 2:
    varValrange = var[irange[0]:rangedim(erange[0], varshape[0]),                               \
      irange[1]:rangedim(erange[1]), varshape[1]]
  elif Nvarshape == 3:
    varValrange = var[irange[0]:rangedim(erange[0], varshape[0]),                               \
      irange[1]:rangedim(erange[1], varshape[1]), irange[2]:rangedim(erange[2], varshape[2])]
  elif Nvarshape == 4:
    varValrange = var[irange[0]:rangedim(erange[0], varshape[0]),                               \
      irange[1]:rangedim(erange[1], varshape[1]), irange[2]:rangedim(erange[2], varshape[2]),      \
      irange[3]:rangedim(erange[3], varshape[3])]
  elif Nvarshape == 5:
    varValrange = var[irange[0]:rangedim(erange[0], varshape[0]),                               \
      irange[1]:rangedim(erange[1], varshape[1]), irange[2]:rangedim(erange[2], varshape[2]),      \
      irange[3]:rangedim(erange[3], varshape[3]), irange[4]:rangedim(erange[4], varshape[4])]
  elif Nvarshape == 6:
    varValrange = var[irange[0]:rangedim(erange[0], varshape[0]),                               \
      irange[1]:rangedim(erange[1], varshape[1]), irange[2]:rangedim(erange[2], varshape[2]),      \
      irange[3]:rangedim(erange[3], varshape[3]), irange[4]:rangedim(erange[4], varshape[4]),      \
      irange[5]:rangedim(erange[5], varshape[5])]

  ncf.close()

  Nshaperange = len(varValrange.shape)
  Noneshape = 0
  for ir in range(Nshaperange):
    if varValrange.shape[ir] == 1:
      Noneshape = Noneshape + 1

  if Noneshape == Nshaperange - 1:
    for i in range(len(varValrange)):
      print '%2s %f' % ( 'NC', varValrange[i] )

def chdimname(values, ncfile, varn):
  """ Changing the name of the dimension
  values = [olddimname]:[newdimname]
    [olddimname]: old name of the dimension
    [newdimname]: new name of the dimension
  ncfile = netCDF file name
  varn = variable name
  """

  if not os.path.isfile(ncfile):
    print errormsg
    print '    chdimname: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')

  olddimname=values.split(':')[0]
  newdimname=values.split(':')[1]

  if not ncf.dimensions.has_key(olddimname):
    print warnmsg
    print '    chdimname: File "' + ncfile + '" does not have dimension "' + olddimname + '" !!!!'
    ncf.close()
    quit()

  if not olddimname == newdimname and not ncf.dimensions.has_key(newdimname):
      newname = ncf.renameDimension(olddimname, newdimname)
      ncf.sync()
      ncf.close()
  else:
      print warnmsg
      print '    chdimname: File "' + ncfile + '" already has dimension name "' + newdimname + '" '
      print '    chdimname: modifying all the variables which use the old dimension'
      filevars = ncf.variables
      for fvarn in filevars:
          if ncf.variables.has_key(fvarn):
              fvar = ncf.variables[fvarn]
              fvardims = fvar.dimensions
              if searchInlist(fvardims, olddimname):
                  print '    variable "' + fvarn + '" uses dimension "' + olddimname + '" '
                  varinf = variable_inf(fvar)

                  newdims = tuple(fvardims)
                  change = {olddimname: newdimname}
# From http://stackoverflow.com/questions/9067043/python-replace-list-values-using-a-tuple
                  newdims = tuple([ change.get(x,x) for x in fvardims ])
                  newvar = ncf.createVariable(fvarn + 'tmpdname', varinf.dtype, newdims, fill_value=varinf.FillValue)
                  varv = fvar[:]
                  newvar[:] = varv

                  for attrn in varinf.attributes:
                      attrv = fvar.getncattr(attrn)
                      newar = newvar.setncattr(attrn, attrv)

                  fvar = ncf.renameVariable(fvarn, fvarn + 'rmdname')
                  ncf.sync()
                  newvar = ncf.renameVariable(fvarn + 'tmpdname', fvarn)
                  ncf.sync()
                  ncf.close()
                  varrm(ncfile, fvarn + 'rmdname')
                  ncf = NetCDFFile(ncfile,'a')

def chvarname(values, ncfile, varn):
  """Changing the name of the variable
  values = new variable name
  ncfile = netCDF file
  varn = name of the variable
  """

  if not os.path.isfile(ncfile):
    print errormsg
    print '    chvarname: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')

  if ncf.dimensions.has_key('plev'):
    # removing pressure level from variable name
    varn = re.sub("\d+", "", varn) 

  if not ncf.variables.has_key(varn):
    print warnmsg
    print '    chvarname: File does not have variable "' + varn + '" !!!!'
    ncf.close()
    quit(-1)

  if not varn == values and not ncf.variables.has_key(values):
    newname = ncf.renameVariable(varn, values)

  ncf.sync()
  ncf.close()

def searchInlist(listname, nameFind):
    """ Function to search a value within a list
    listname = list
    nameFind = value to find
    >>> searInlist(['1', '2', '3', '5'], '5')
    True
    """
    for x in listname:
      if x == nameFind:
        return True
    return False

def set_attribute(ncvar, attrname, attrvalue):
    """ Sets a value of an attribute of a netCDF variable. Removes previous attribute value if exists
    ncvar = object netcdf variable
    attrname = name of the attribute
    attrvalue = value of the attribute
    """
    import numpy as np
    from netCDF4 import Dataset as NetCDFFile

    attvar = ncvar.ncattrs()
    if searchInlist(attvar, attrname):
        attr = ncvar.delncattr(attrname)

    attr = ncvar.setncattr(attrname, attrvalue)

    return ncvar

def set_attributek(ncvar, attrname, attrval, attrkind):
    """ Sets a value of an attribute of a netCDF variable with a kind. Removes previous attribute value if exists
    ncvar = object netcdf variable
    attrname = name of the attribute
    attrvalue = value of the attribute
    atrtrkind = kind of attribute: 'S', string ('!' as spaces); 'I', integer; 'R', real; 'D', double
    """
    import numpy as np
    from netCDF4 import Dataset as NetCDFFile

    if attrkind == 'S':
        attrvalue = str(attrval.replace('!', ' '))
    elif attrkind == 'I':
        attrvalue = int(attrval)
    elif attrkind == 'R':
        attrvalue = float(attrval)
    elif attrkind == 'D':
        attrvalue = np.float64(attrval)
    else:
        print errormsg
        print '    set_attributek: "' + attrkind + '" kind of attribute is not ready!'
        quit(-1)

    attvar = ncvar.ncattrs()
    if searchInlist(attvar, attrname):
        attr = ncvar.delncattr(attrname)
    attr = ncvar.setncattr(attrname, attrvalue)

    return ncvar

def gaddattr(values, ncfile):
  """ Add a global attribute to a netCDF. Removes previous attribute if it exist
  values = [attrname]|[attrvalue ('!' as spaces)]
  ncfile = netCDF file
  """
  if not os.path.isfile(ncfile):
    print errormsg
    print '    gaddattr: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')

  attrvals=values.split('|')
  attrn=attrvals[0]
  attrv=attrvals[1].replace('!', ' ')

  ncf = set_attribute(ncf, attrn, attrv)

  ncf.sync()
  ncf.close()

def gaddattrk(values, ncfile):
  """ Add a global attribute to a netCDF caring about the type. Removes previous attribute if it exist
  values = [attrname]|[attrvalue]|[attrk]
    attrname = name of the attribute
    attrvalue = value of the attribute
    attrk = 'S', string ('!' as spaces); 'I', integer; 'R', real; 'D', double
  ncfile = netCDF file
  """
  if not os.path.isfile(ncfile):
    print errormsg
    print '    gaddattrk: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')

  attrvals=values.split('|')
  attrn=attrvals[0]
  attrv0=attrvals[1]
  attrk=attrvals[2]

  if attrk == 'S':
    attrv = str(attrv0.replace('!', ' '))
  elif attrk == 'I':
    attrv = int(attrv0)
  elif attrk == 'R':
    attrv = float(attrv0)
  elif attrk == 'D':
    attrv = np.float64(attrv0)
  else:
    print errormsg
    print '    gaddattrk: "' + attrk + '" kind of attribute is not ready!'
    ncf.close()
    quit(-1)

  ncf = set_attribute(ncf, attrn, attrv)

  ncf.sync()
  ncf.close()

def varaddattr(values, ncfile, varn):
  """ Add an attribute to a variable. Removes previous attribute if it exists
  values = [attrname]|[attrvalue('!' as spaces)]
  ncfile = netCDF file name
  varn = name of the variable
  """
  if not os.path.isfile(ncfile):
    print errormsg
    print '    varaddattr: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')

  if ncf.dimensions.has_key('plev'):
    # removing pressure level from variable name
    varn = re.sub("\d+", "", varn) 

  if not ncf.variables.has_key(varn):
    print errormsg
    print '    varaddattr: File "' + ncfile + '" does not have variable "' + varn + '" !!!!'
    print errormsg
    ncf.close()
    quit(-1)

  attrvals=values.split('|')
  attrn=attrvals[0]
  attrv=attrvals[1].replace('!', ' ')

  var = ncf.variables[varn]
  var = set_attribute(var, attrn, attrv)

  ncf.sync()
  ncf.close()

def varaddattrk(values, ncfile, varn):
  """ Add an attribute to a variable caring about the type
  values = [attrname]|[attrvalue]|[attrk]
    attrname = name of the attribute
    attrvalue = value of the attribute
    attrk = 'S', string ('!' as spaces); 'I', integer; 'R', real; 'D', double
  ncfile = netCDF file
  varn = variable name
  """
  if not os.path.isfile(ncfile):
    print errormsg
    print '    varaddattrk: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')

  if ncf.dimensions.has_key('plev'):
    # removing pressure level from variable name
    varn = re.sub("\d+", "", varn) 

  if not ncf.variables.has_key(varn):
    print errormsg
    print '    varaddattrk: File "' + ncfile + '"does not have variable ' + varn + ' !!!!'
    print errormsg
    ncf.close()
    quit(-1)

  attrvals=values.split('|')
  attrn=attrvals[0]
  attrv0=attrvals[1]
  attrk=attrvals[2]

  var = ncf.variables[varn]
  if attrk == 'S':
    attrv = str(attrv0.replace('!', ' '))
  elif attrk == 'I':
    attrv = int(attrv0)
  elif attrk == 'R':
    attrv = float(attrv0)
  elif attrk == 'D':
    attrv = np.float64(attrv0)
  else:
    print errormsg
    print '    varaddattrk: "' + attrk + '" kind of attribute is not ready!'
    ncf.close()
    quit(-1)

  var = set_attribute(var, attrn, attrv)

  ncf.sync()
  ncf.close()

def varrmattr(values, ncfile, varn):
  """ Removing an attribute from a variable
  values = attribute name
  ncfile = netCDF file name
  varn = variable name
  """
  if not os.path.isfile(ncfile):
    print errormsg
    print '    varrmattr: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')

  if ncf.dimensions.has_key('plev'):
    # removing pressure level from variable name
    varn = re.sub("\d+", "", varn) 

  if not ncf.variables.has_key(varn):
    print errormsg
    print '    varrmattr: File "' + ncfile + '" does not have variable "' + varn + '" !!!!'
    print errormsg
    ncf.close()
    quit(-1)

  var = ncf.variables[varn]

  attvar = var.ncattrs()
  if searchInlist(attvar, values):
      attr = var.delncattr(values)
  else:
      print warnmsg
      print '    varrmattr: "' + varn + '" does not have attribute: ' + values

  ncf.sync()
  ncf.close()

def grmattr(values, ncfile):
  """ Removing a global attribute
  values = attribute name
  ncfile = netCDF file
  """
  if not os.path.isfile(ncfile):
    print errormsg
    print '    grmattr: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')

  attvar = ncf.ncattrs()
  if searchInlist(attvar, values):
      attr = ncf.delncattr(values)
  else:
      print warnmsg
      print '  grmattr: "' + ncfile + '" does not have attribute: ' + values

  ncf.sync()
  ncf.close()

def fvaradd(values, ncfile):
  """ Adding variable (and all its attributes and dimensions) from a reference file to a file
  values = [netCDFref]:[varnref]
    netCDFref = netCDF file name as reference for the variable to add
    varnref = name of the variable from [netCDFref] to be added
  ncfile = netCDF file name
  """

  refnc = values.split(':')[0]
  refvar = values.split(':')[1]

  if not os.path.isfile(ncfile):
    print errormsg
    print '    fvaradd: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  if not os.path.isfile(refnc):
    print errormsg
    print '    fvaradd: Reference file "' + refnc + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')
  ncref = NetCDFFile(refnc,'r')

  refvars = ncref.variables
  if searchInlist(refvars, refvar):
      refvarv = ncref.variables[refvar]
  else:
      print errormsg
      print '    fvaradd: File "' + refnc + '" does not have variable: ' + refvar
      ncf.close()
      ncref.close()
      quit(-1)

  vardims = refvarv.dimensions
  vartype = refvarv.dtype
  varattr = refvarv.ncattrs()

# Checking dimensions
##
  newdims = ncf.dimensions
  for rdim in vardims:
      if not searchInlist(newdims, rdim):
          print '      fvaradd: Adding dimension ' + rdim
          ncf.close()
          ncref.close()
          fdimadd(refnc + ':' + rdim, ncfile)
          ncf = NetCDFFile(ncfile,'a')
          ncref = NetCDFFile(refnc,'r')

# Checking fill value
## 
  if searchInlist(varattr, '_FillValue'):
      varfil = refvarv._FillValue
  else:
      varfil = False

  print '      fvaradd: Adding refvar:', refvar, 'shape: ', refvarv.shape
  var = ncf.createVariable(refvar, vartype, vardims, fill_value=varfil)

  if not len(refvarv.shape) == 0:
    var[:] = refvarv[:]

  newvar = ncf.variables[refvar]
  for attr in varattr:
      newvarattrs = newvar.ncattrs()
      attrv = refvarv.getncattr(attr)
      if not searchInlist(newvarattrs, attr):     
          newvar.setncattr(attr, attrv)

  ncf.sync()
  ncf.close()
  ncref.close()

def fdimadd(values, ncfile):
  """ Adding dimension from another reference file 
  values = [refnc]:[refdim]
    refnc = netCDF file name as reference
    refdim = name of the dimension to be added
  ncfile = netCDF file name
  """

  refnc = values.split(':')[0]
  refdim = values.split(':')[1]

  if not os.path.isfile(ncfile):
    print errormsg
    print '    fdimadd: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  if not os.path.isfile(refnc):
    print errormsg
    print '    fdimadd: Reference file "' + refnc + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')
  ncref = NetCDFFile(refnc,'r')

  refdims = ncref.dimensions
  if searchInlist(refdims, refdim):
      refdimv = ncref.dimensions[refdim]
  else:
      print errormsg
      print '    fdimadd: File "' + refnc + '" does not have dimension: "' + refdim + '"'
      ncf.close()
      ncref.close()
      quit(-1)

  if refdimv.isunlimited():
      print '      fdimadd: Unlimited dimension '
      dimsize = None
  else:
      dimsize = len(refdimv)

  print '      fdimadd: Adding refdim:', refdim, 'size:', dimsize
  dim = ncf.createDimension(refdim, dimsize)
  
  ncf.sync()
  ncf.close()
  ncref.close()

def fattradd(var, values, ncfile):
  """ Adding attributes from a reference file
  var = variable to which has to be added the attribute
  values = [refnc]:[refvar]
    refnc = netCDF file name as reference
    refvar = variable from the reference file
  ncfile = netCDF file name
  """

  refnc = values.split(':')[0]
  refvar = values.split(':')[1]

  if not os.path.isfile(ncfile):
    print errormsg
    print '    fattradd: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  if not os.path.isfile(refnc):
    print errormsg
    print '    fattradd: Reference file "' + refnc + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')
  ncref = NetCDFFile(refnc,'r')

  vars = ncf.variables
  if searchInlist(vars, var):
      varv = ncf.variables[var]
  else:
      print '  fattradd: File "' + ncfile + '" does not have variable: "' + var + '"'
      ncf.close()
      ncref.close()
      quit(-1)

  refvars = ncref.variables
  if searchInlist(refvars, refvar):
      refvarv = ncref.variables[refvar]
  else:
      print '    fattradd: File "' + refnc + '" does not have variable: "' + refvar + '"'
      ncf.close()
      ncref.close()
      quit(-1)

  refvarattrs = refvarv.ncattrs()
  Nattrs = len(refvarattrs)
  print '      fattradd: Adding ', Nattrs,' atributes from:', refvar

  for attr in refvarattrs:
      attrv = refvarv.getncattr(attr)
      atvar = set_attribute(varv, attr, attrv)
  
  ncf.sync()
  ncf.close()
  ncref.close()

def fgaddattr(values, ncfile):
  """ Adding global attributes from a reference file
  values = netCDF file name as reference
  ncfile = netCDF file name
  """

  refnc = values

  if not os.path.isfile(ncfile):
    print errormsg
    print '    fgaddattr: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  if not os.path.isfile(refnc):
    print errormsg
    print '    fgaddattr: Reference file "' + refnc + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')
  ncref = NetCDFFile(refnc,'r')

  refgattrs = ncref.ncattrs()
  Nattrs = len(refgattrs)
  print '      fgaddattr: Adding ', Nattrs,' global atributes'

  for attr in refgattrs:
      attrv = ncref.getncattr(attr)
      atvar = set_attribute(ncf, attr, attrv)
  
  ncf.sync()
  ncf.close()
  ncref.close()

def varrm(ncfile, var):
  """ Removing a variable from a file
  ncfile = netCDF file name
  var = variable name to remove
  """
  import shutil as shu

  if not os.path.isfile(ncfile):
    print errormsg
    print '    varrm: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')
  ncvars = ncf.variables
  ncf.close()

  if not searchInlist(ncvars, var):
      print '    varrm: File "' + ncfile + '" does not have variable: "' + var + '"'
      ncf.close()
      quit(-1)

  tmpncf = NetCDFFile('tmp_py.nc' , 'w')
  gtmpattr = set_attribute(tmpncf, 'copy', 'temporal')
  tmpncf.sync()
  tmpncf.close()

  for varn in ncvars:
      if not varn == var:
           fvaradd(ncfile + ':' + varn, 'tmp_py.nc')

  fgaddattr(ncfile, 'tmp_py.nc')
  shu.copyfile('tmp_py.nc', ncfile)
  os.remove('tmp_py.nc')
 
def ivars(ncfile):
  """Give all the variable names of a file
  """
  if not os.path.isfile(ncfile):
    print errormsg
    print '    ivars: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'a')
  ncvars = ncf.variables
  allvars=''
  for var in ncvars:
      print var
      allvars=allvars + ':' + var

  print '  # allvars= ' + allvars
  ncf.close()

def igattrs(ncfile):
  """Give all the global attributes of a file
  """
  if not os.path.isfile(ncfile):
    print errormsg
    print '    igattrs: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit()    

  ncf = NetCDFFile(ncfile,'r')
  ncattrs = ncf.ncattrs()
  allattrs=''
  for attr in ncattrs:
      attrv=ncf.getncattr(attr)
      if type(attrv) == type(str('1')):
          attrk = 'S'
      elif type(attrv) == type(unicode('1')):
          attrk = 'S'
      elif type(attrv) == type(int(1)):
          attrk = 'I'
      elif type(attrv) == type(np.int(1)):
          attrk = 'I'
      elif type(attrv) == type(np.int32(1)):
          attrk = 'I'
      elif type(attrv) == type(float(1.)):
          attrk = 'R'
      elif type(attrv) == type(np.float32(1.)):
          attrk = 'R'
      elif type(attrv) == type(np.float64(1.)):
          attrk = 'D'
      else:
          print errormsg
          print '    igattr: Reading attribute "', type(attrv), '" not ready! value:', attrv
          ncf.close()
          quit(-1)
      print attr, '|',  attrv, '|', attrk
      allattrs=allattrs + ':' + attr + '|' + str(attrv) + '|' + attrk

  print '####### ###### ##### #### ### ## #'
  print '# allgattrs= ' + allattrs
  ncf.close()

def isgattrs(values, ncfile):
  """Give a single global attribute of a file and its type
  values = attribute name
  ncfile = netCDF file name
  output:
    attribute name, '|',  attribute value, '|', attribute kind ('S', string '!' as spaces; 'I', integer; 'R', real; 'D', double )
  """

  if not os.path.isfile(ncfile):
    print errormsg
    print '    isgattrs: File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'r')
  ncattrs = ncf.ncattrs()

  if not searchInlist(ncattrs, values):
      print '    isgattrs: File "' + ncfile + '" does not have global attribute: "' + values + '" '
      ncf.close()
      quit(-1)

  attrv=ncf.getncattr(values)
  if type(attrv) == type(str('1')):
      attrk = 'S'
      attrv = attrv.replace(' ','!')
  elif type(attrv) == type(unicode('1')):
      attrk = 'S'
      attrv = attrv.replace(' ','!')
  elif type(attrv) == type(int(1)):
      attrk = 'I'
  elif type(attrv) == type(np.int(1)):
      attrk = 'I'
  elif type(attrv) == type(np.int32(1)):
      attrk = 'I'
  elif type(attrv) == type(float(1.)):
      attrk = 'R'
  elif type(attrv) == type(np.float32(1.)):
      attrk = 'R'
  elif type(attrv) == type(np.float64(1.)):
      attrk = 'D'
  else:
      print errormsg
      print '    isgattr: Reading attribute "', type(attrv), '" not ready! value:', attrv
      ncf.close()
      quit(-1)
  print values, '|',  attrv, '|', attrk

  ncf.close()

def ivattrs(ncfile, varn):
  """Give all the attributes of a variable and its type
  ncfile = netCDF file name
  var = variable name
  output:
    attribute name, '|',  attribute value, '|', attribute kind ('S', string '!' as spaces; 'I', integer; 'R', real; 'D', double )
  """

  if not os.path.isfile(ncfile):
    print errormsg
    print '    ivattrs:File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'r')
  ncvars = ncf.variables

  if not searchInlist(ncvars, varn):
      print '    ivattrs:"' + ncfile + '" does not have variable: "' + varn + '" '
      ncf.close()
      quit(-1)
  varval = ncf.variables[varn]

  ncattrs = varval.ncattrs()
  allattrs=''
  for attr in ncattrs:
      attrv=varval.getncattr(attr)
      if type(attrv) == type(str('1')):
          attrk = 'S'
          attrv = attrv.replace(' ','!')
      elif type(attrv) == type(unicode('1')):
          attrk = 'S'
          attrv = attrv.replace(' ','!')
      elif type(attrv) == type(int(1)):
          attrk = 'I'
      elif type(attrv) == type(np.int(1)):
          attrk = 'I'
      elif type(attrv) == type(np.int32(1)):
          attrk = 'I'
      elif type(attrv) == type(float(1.)):
          attrk = 'R'
      elif type(attrv) == type(np.float32(1.)):
          attrk = 'R'
      elif type(attrv) == type(np.float64(1.)):
          attrk = 'D'
      else:
          print errormsg
          print '    ivattrs: Reading attribute "', type(attrv), '" not ready! value:', attrv
          ncf.close()
          quit(-1)
      print attr, '|',  attrv, '|', attrk
      allattrs=allattrs + ':' + attr + '|' + str(attrv) + '|' + attrk

  print '####### ###### ##### #### ### ## #'
  print '# allvattrs= ' + allattrs
  ncf.close()

def isvattrs(values, ncfile, varn):
  """Give a single attribute of a variable
  values = attribute name
  ncfile = netCDF file name
  varn = variable name
  """

  if not os.path.isfile(ncfile):
    print errormsg
    print '    isvattrs:File "' + ncfile + '" does not exist !!'
    print errormsg
    quit(-1)    

  ncf = NetCDFFile(ncfile,'r')
  ncvars = ncf.variables

  if not searchInlist(ncvars, varn):
      print errormsg
      print '    isvattrs:"' + ncfile + '" does not have variable: "' + varn + '" '
      ncf.close()
      quit(-1)
  varval = ncf.variables[varn]

  ncattrs = varval.ncattrs()
  if not searchInlist(ncattrs, values):
      print errormsg
      print '    isvattrs:' + ncfile + ' does not have global attribute: "' + values + '" '
      ncf.close()
      quit(-1)

  attrv=varval.getncattr(values)
  if type(attrv) == type(str('1')):
      attrk = 'S'
  elif type(attrv) == type(unicode('1')):
      attrk = 'S'
  elif type(attrv) == type(int(1)):
      attrk = 'I'
  elif type(attrv) == type(np.int(1)):
      attrk = 'I'
  elif type(attrv) == type(np.int32(1)):
      attrk = 'I'
  elif type(attrv) == type(float(1.)):
      attrk = 'R'
  elif type(attrv) == type(np.float(1.)):
      attrk = 'R'
  elif type(attrv) == type(np.float32(1.)):
      attrk = 'R'
  elif type(attrv) == type(np.float64(1.)):
      attrk = 'D'
  else:
      print errormsg
      print '    isvattr: Reading attribute "', type(attrv), '" not ready! value:', attrv
      ncf.close()
      quit(-1)
  print values, '|',  attrv, '|', attrk

  ncf.close()

def grattr(values, ncfile):
  """ Function to read a global atribute
  values = attribute name
  ncfile = netCDF file name
  """
  ncf = NetCDFFile(ncfile,'r')

  glob_attrs = ncf.ncattrs()

  attrPos = searchInlist(glob_attrs, values)

  if not attrPos:
    print errormsg
    print '    grattr: File "' + ncfile + '" does not have attribute "' + values + '" !!!!'
    print errormsg
    ncf.close()
    quit(-1)

  print ncf.getncattr(values)

  ncf.close()

def vrattr(values, ncfile, varn):
  """ Function to remove an atribute from a variable
  values = attribute name
  ncfile = netCDF file name
  varn = variable name
  """
  ncf = NetCDFFile(ncfile,'r')

  if not ncf.variables.has_key(varn):
    print errormsg
    print '    vrattr: File "' + ncfile + '" does not have variable "' + varn + '" !!!!'
    print errormsg
    ncf.close()
    quit(-1)

  var = ncf.variables[varn]
  var_attrs = var.ncattrs()

  attrPos = searchInlist(var_attrs, values)

  if not attrPos:
    print errormsg
    print '   vrattr: Variable "' + varn + '" does not have attribute "' + values + '" !!!!'
    print errormsg
    ncf.close()
    quit(-1)

  print var.getncattr(values)

  ncf.close()

def datetimeStr_datetime(StringDT):
  """ Function to transform a string date ([YYYY]-[MM]-[DD]_[HH]:[MI]:[SS] format ) to a date object
  >>> datetimeStr_datetime('1976-02-17_00:00:00')
  1976-02-17 00:00:00
  """
  import datetime as dt

  dateDT = StringDT.split('_')
  dateD = dateDT[0].split('-')
  timeT = dateDT[1].split(':')

  if int(dateD[0]) == 0:
    print warnmsg
    print '    datetimeStr_datetime: 0 reference year!! changing to 1'
    dateD[0] = 1 

  if len(timeT) == 3:
    newdatetime = dt.datetime(int(dateD[0]), int(dateD[1]), int(dateD[2]), int(timeT[0]), int(timeT[1]), int(timeT[2]))
  else:
    newdatetime = dt.datetime(int(dateD[0]), int(dateD[1]), int(dateD[2]), int(timeT[0]), int(timeT[1]), 0)

  return newdatetime

def dateStr_date(StringDate):
  """ Function to transform a string date ([YYYY]-[MM]-[DD] format) to a date object
  >>> dateStr_date('1976-02-17')
  1976-02-17
  """
  import datetime as dt

  dateD = StringDate.split('-')
  if int(dateD[0]) == 0:
    print warnmsg
    print '    dateStr_date: 0 reference year!! changing to 1'
    dateD[0] = 1
  newdate = dt.date(int(dateD[0]), int(dateD[1]), int(dateD[2]))
  return newdate

def timeStr_time(StringDate):
  """ Function to transform a string date ([HH]:[MI]:[SS] format) to a time object
  >>> datetimeStr_datetime('04:32:54')
  04:32:54
  """
  import datetime as dt

  timeT = StringDate.split(':')
  if len(timeT) == 3:
    newtime = dt.time(int(timeT[0]), int(timeT[1]), int(timeT[2]))
  else:
    newtime = dt.time(int(timeT[0]), int(timeT[1]), 0)

  return newtime

def time_information(ncfu, tname):
    """ Function to provide information about variable time
    ncfu = netCDF unit name
    tname = name of the variable time in [ncfu]
    """
    times = ncfu.variables[tname]
    timeinf = []

    attvar = times.ncattrs()
    if not searchInlist(attvar, 'units'):
        print errormsg
        print '    "time" does not have attribute: "units"'
        quit(-1)
    else:
        units = times.getncattr('units')
  
    txtunits = units.split(' ')
    tunits = txtunits[0]
    Srefdate = txtunits[len(txtunits) - 1]
# Does reference date contain a time value [YYYY]-[MM]-[DD] [HH]:[MI]:[SS]
##
    timeval = Srefdate.find(':')

    if not timeval == -1:
#        print '  refdate with time!'
        refdate = datetimeStr_datetime(txtunits[len(txtunits) - 2] + '_' + Srefdate)
    else:
        refdate = dateStr_date(Srefdate)

    timeinf.append(tunits)
    timeinf.append(Srefdate)
    timeinf.append(refdate)

    return timeinf

def CFtimes_datetime(ncfile, tname):
    """ Provide date/time array from a file with a series of netCDF CF-compilant time variable
    ncfile = netCDF file name
    tname = name of the variable time in [ncfile]
    output:
      array(dimt, 0) = year
      array(dimt, 1) = month
      array(dimt, 2) = day
      array(dimt, 3) = hour
      array(dimt, 4) = minute
      array(dimt, 5) = second
    """
    import datetime as dt

    times = ncfile.variables[tname]
    timevals = times[:]

    attvar = times.ncattrs()
    if not searchInlist(attvar, 'units'):
        print errormsg
        print '    CFtimes_datetime: "time" does not have attribute: "units"'
        quit(-1)
    else:
        units = times.getncattr('units')
  
    txtunits = units.split(' ')
    tunits = txtunits[0]
    Srefdate = txtunits[len(txtunits) - 1]
# Does reference date contain a time value [YYYY]-[MM]-[DD] [HH]:[MI]:[SS]
##
    timeval = Srefdate.find(':')

    if not timeval == -1:
#        print '  refdate with time!'
        refdate = datetimeStr_datetime(txtunits[len(txtunits) - 2] + '_' + Srefdate)
    else:
        refdate = dateStr_date(Srefdate)

    dimt = len(timevals)
    realdates = np.zeros((dimt, 6), dtype=int)
    print realdates.shape

## Not in timedelta
#    if tunits == 'years':
#        for it in range(dimt):
#            realdate = refdate + dt.timedelta(years=float(times[it]))
#            realdates[it] = int(realdate.year)
#    elif tunits == 'months':
#        for it in range(dimt):
#            realdate = refdate + dt.timedelta(months=float(times[it]))
#            realdates[it] = int(realdate.year)
    if tunits == 'weeks':
        for it in range(dimt):
            realdate = refdate + dt.timedelta(weeks=float(times[it]))
            realdates[it,0] = int(realdate.year)
            realdates[it,1] = int(realdate.month)
            realdates[it,2] = int(realdate.day)
            realdates[it,3] = int(realdate.hour)
            realdates[it,4] = int(realdate.second)
            realdates[it,5] = int(realdate.minute)
    elif tunits == 'days':
        for it in range(dimt):
            realdate = refdate + dt.timedelta(days=float(times[it]))
            realdates[it,0] = int(realdate.year)
            realdates[it,1] = int(realdate.month)
            realdates[it,2] = int(realdate.day)
            realdates[it,3] = int(realdate.hour)
            realdates[it,4] = int(realdate.second)
            realdates[it,5] = int(realdate.minute)
    elif tunits == 'hours':
       for it in range(dimt):
            realdate = refdate + dt.timedelta(hours=float(times[it]))
            realdates[it,0] = int(realdate.year)
            realdates[it,1] = int(realdate.month)
            realdates[it,2] = int(realdate.day)
            realdates[it,3] = int(realdate.hour)
            realdates[it,4] = int(realdate.second)
            realdates[it,5] = int(realdate.minute)
    elif tunits == 'minutes':
       for it in range(dimt):
            realdate = refdate + dt.timedelta(minutes=float(times[it]))
            realdates[it,0] = int(realdate.year)
            realdates[it,1] = int(realdate.month)
            realdates[it,2] = int(realdate.day)
            realdates[it,3] = int(realdate.hour)
            realdates[it,4] = int(realdate.second)
            realdates[it,5] = int(realdate.minute)
    elif tunits == 'seconds':
       for it in range(dimt):
            realdate = refdate + dt.timedelta(seconds=float(times[it]))
            realdates[it,0] = int(realdate.year)
            realdates[it,1] = int(realdate.month)
            realdates[it,2] = int(realdate.day)
            realdates[it,3] = int(realdate.hour)
            realdates[it,4] = int(realdate.second)
            realdates[it,5] = int(realdate.minute)
    elif tunits == 'milliseconds':
       for it in range(dimt):
            realdate = refdate + dt.timedelta(milliseconds=float(times[it]))
            realdates[it,0] = int(realdate.year)
            realdates[it,1] = int(realdate.month)
            realdates[it,2] = int(realdate.day)
            realdates[it,3] = int(realdate.hour)
            realdates[it,4] = int(realdate.second)
            realdates[it,5] = int(realdate.minute)
    else:
          print errormsg
          print '    CFtimes_datetime: time units "' + tunits + '" not ready!!!!'
          quit(-1)

    return realdates

class variable_inf(object):
    """ Class to provide the information from a given variable
    var = object netCDF variable
    self.name: name of the variable
    self.dtype: type of the variable
    self.attributes: list with the name of attributes
    self.FillValue: value of the missing value
    self.dims: dimensions of the variable
    self.Ndims: number of dimensions
    self.dimx: length of dimension along x-axis
    self.dimy: length of dimension along y-axis
    self.sname: standard name
    self.lname: long name
    self.corr: attribute 'coordinates'
    self.units: units of the variable
    """

    def __init__(self, var):

        if var is None:
            self.name = None
            self.dims = None
            self.Ndims = None
            self.dimx = None
            self.dimy = None
            self.sname = None
            self.lname = None
            self.corr = None
            self.units = None
            self.FillValue = None
            self.dtype = None
            self.attributes = None
        else:
            self.name = var._name
            self.dtype = var.dtype
            self.attributes = var.ncattrs()
            self.dims = var.shape

            if searchInlist(self.attributes, 'standard_name'):
                self.sname = var.getncattr('standard_name')
            else:
                print '    variable_inf.classpy: variable "' + self.name + '" does not have attribute "standard_name"'
                self.sname = None

            if searchInlist(self.attributes, 'long_name'):
                self.lname = var.getncattr('long_name')
            else:
                print '    variable_inf.classpy: variable "' + self.name + '" does not have attribute "long_name"'
                self.lname = None

            if searchInlist(self.attributes, 'coordinates'):
                self.coor = var.getncattr('coordinates')
            else:
                print '    variable_inf.classpy: variable "' + self.name + '" does not have attribute "coordinates"'
                self.coor = None

            if searchInlist(self.attributes, 'units'):
                self.units = var.getncattr('units')
            else:
                print '    variable_inf.classpy: variable "' + self.name + '" does not have attribute "units"'
                self.units = None

            if searchInlist(self.attributes, '_FillValue'):
                self.FillValue = var.getncattr('_FillValue')
            else:
                print '    variable_inf.classpy: variable "' + self.name + '" does not have attribute "_FillValue"'
                self.FillValue = None
             
            self.Ndims = len(self.dims)
            if self.Ndims == 1:
                self.dimx=self.dims[0]
            if self.Ndims == 2:
                self.dimy=self.dims[0]
                self.dimx=self.dims[1]
            if self.Ndims == 3:
                self.dimy=self.dims[1]
                self.dimx=self.dims[2]
            if self.Ndims == 4:
                self.dimy=self.dims[2]
                self.dimx=self.dims[3]

def subyrs(values, ncfile, varn):
  """ Function to retrieve a series of years from a file
  values = 
    [year1]:[[year2]:...[yearn]] values for years [year1], [year2], ... [yearn]
    [yearI]-[yearE] values for the period between [yearI] and [yearN]
  ncfile = netCDF file name
  varn = variable name
  """
  import datetime as dt
  import calendar as cal
  ofile = 'subsetyrs.nc'

  ncf = NetCDFFile(ncfile,'r')

  if not ncf.variables.has_key(varn):
    print errormsg
    print '   subyrs: File does not have variable "' + varn + '" !!!!'
    print errormsg
    ncf.close()
    quit(-1)

  times = ncf.variables['time']
  timevals = times[:]
  realdates = CFtimes_datetime(ncf, 'time')
  dimt = len(timevals)

# Checking years
##
  isper = values.find('-')

  desiredvalues = np.array(dimt, dtype=bool)
  desiredvaluesi = np.array(dimt, dtype=bool)
  desiredvaluese = np.array(dimt, dtype=bool)

  if not isper == -1:
# Years list given as a period [YYYYi]-[YYYYf]
    print '      subyrs: There is a period of years "' + values + '"'
    iyr = int(values.split('-')[0])
    eyr = int(values.split('-')[1])
    desiredvaluesi = np.array(realdates[:,0] >= iyr)
    desiredvaluese = np.array(realdates[:,0] <= eyr)

    desiredvalues = desiredvaluesi*desiredvaluese

  else:
    yrs = values.split(':')
    nyr = 1
    for iyr in range(len(yrs)):
      desiredvaluesi = np.array(realdates[:,0] == int(yrs[iyr]))
      if nyr == 1:
        desiredvalues = desiredvaluesi
      else:
        desiredvalues = desiredvaluesi+desiredvalues
      nyr = nyr + 1

  Ndesiredvalues = len(realdates[desiredvalues])

  print '      subyrs: N values: ', Ndesiredvalues
  if Ndesiredvalues == 0:
    print errormsg
    print '    subyrs: No values found for "' + values + '"!!'
    ncf.close()
    quit(-1)

# Variable values (assuming time as first dimension)
##
  var = ncf.variables[varn]
#  varvals = var[:]
  vardims = var.shape
  varshape = len(vardims)
  print '      subyrs: Shape of data: ',varshape, ':' , vardims

  if varshape == 1:
    vardesiredvalues = np.arange(Ndesiredvalues)
    vardesiredvalues = var[desiredvalues]
  elif varshape == 2:
    vardesiredvalues = np.arange(Ndesiredvalues*vardims[1]).reshape(Ndesiredvalues,vardims[1])
    vardesiredvalues = var[desiredvalues, :]
  elif varshape == 3:
    vardesiredvalues = np.arange(Ndesiredvalues*vardims[1]*vardims[2]).reshape(Ndesiredvalues,                        \
      vardims[1],vardims[2])
    vardesiredvalues = var[desiredvalues, :, :]
  elif varshape == 4:
    vardesiredvalues = np.arange(Ndesiredvalues*vardims[1]*vardims[2]*vardims[3]).reshape(Ndesiredvalues,             \
      vardims[1],vardims[2],vardims[3])
    vardesiredvalues = var[desiredvalues, :, :, :]
  elif varshape == 5:
    vardesiredvalues = np.arange(Ndesiredvalues*vardims[1]*vardims[2]*vardims[3]*vardims[4]).reshape(Ndesiredvalues,  \
      vardims[1],vardims[2],vardims[3],vardims[4])
    vardesiredvalues = var[desiredvalues, :, :, :, :]
  elif varshape == 6:
    vardesiredvalues = np.arange(Ndesiredvalues*vardims[1]*vardims[2]*vardims[3]*vardims[4]*                          \
      vardims[5]).reshape(Ndesiredvalues,vardims[1],vardims[2],vardims[3],vardims[4],vardims[5])
    vardesiredvalues = var[desiredvalues, :, :, :, :, :]
  else:
     print errormsg
     print '    subyrs: ', varshape, ' shape of matrix not prepared !'
     ncf.close()
     quit(-1)

#  print '  shape of desired values: ', vardesiredvalues.shape
# Creation of file
##

  ncfo = NetCDFFile( ofile, 'w')

  vardims = var.dimensions
  vartype = var.dtype
  varattr = var.ncattrs()

# Checking dimensions
##
  newdims = ncfo.dimensions
  for rdim in vardims:
      if not searchInlist(newdims, rdim):
          if not rdim == 'time':
              print '      subyrs: Adding dimension ' + rdim
              ncfo.sync()
              ncf.close()
              ncfo.close()
              fdimadd(ncfile + ':' + rdim, ofile)
              ncf = NetCDFFile(ncfile,'r')
              ncfo = NetCDFFile(ofile,'a')
          else:
              ncfo.createDimension('time', None)

# Checking fill value
## 
  if searchInlist(varattr, '_FillValue'):
      varfil = var._FillValue
  else:
      varfil = False

  newvar = ncfo.createVariable(varn, vartype, vardims, fill_value=varfil)

  if not varshape == 0:
    newvar[:] = vardesiredvalues

  newvar = ncfo.variables[varn]
  for attr in varattr:
      newvarattrs = newvar.ncattrs()
      attrv = var.getncattr(attr)
      if not searchInlist(newvarattrs, attr):     
          newvar.setncattr(attr, attrv)

  vardims = times.dimensions
  vartype = times.dtype
  varattr = times.ncattrs()

  newvar = ncfo.createVariable('time', vartype, vardims, fill_value=varfil)
  newvar = ncfo.variables['time']
  newvar[:] = timevals[desiredvalues]

  ncf.close()
  ncfo.sync()
  ncfo.close()
  fattradd('time', ncfile + ':time', ofile)
  fvaradd(ncfile + ':lon', ofile)
  fvaradd(ncfile + ':lat', ofile)

  ncfo = NetCDFFile(ofile,'a')
  newvar = ncfo.variables['time']
  newvarattrs = newvar.ncattrs()
  if searchInlist(newvarattrs, 'bounds'):
      if newvar.getncattr('bounds') == 'time_bnds':
          ncf = NetCDFFile(ncfile,'r')
          tbnds = ncf.variables['time_bnds']
          vardims = tbnds.dimensions
          vartype = tbnds.dtype
          varattr = tbnds.ncattrs()
          ncfo.createDimension('bnds', 2)
          newvar = ncfo.createVariable('time_bnds', vartype, vardims, fill_value=varfil)
          newvar[:] = tbnds[desiredvalues,:]

          ncf.close()
          ncfo.sync()
          ncfo.close()
          fattradd('time_bnds', ncfile + ':time_bnds', ofile)
      else:
          ncfo.close()
  else:
      ncfo.close()

  fgaddattr(ncfile, ofile)

  print '      subyrs: File "' + ofile + '" with a subset of ' + values + ' has been created'

def submns(values, ncfile, varn):
  """ Function to retrieve a series of months from a file
  values = 
    [mon1]:[[mion2]:...[monn]] values for months [mon1], [mon2], ... [monn]
    [monI]-[monE] values for the period between [monI] and [monN]
  ncfile = netCDF file name
  varn = variable name
  """
  import datetime as dt
  import calendar as cal

  ofile = 'subsetmns.nc'

  ncf = NetCDFFile(ncfile,'r')

  if not ncf.variables.has_key(varn):
    print errormsg
    print '    submns: File "' + ncfile + '" does not have variable "' + varn + '" !!!!'
    print errormsg
    ncf.close()
    quit(-1)

  times = ncf.variables['time']
  timevals = times[:]
  realdates = CFtimes_datetime(ncf, 'time')
  dimt = len(timevals)

# Checking months
##
  isper =values.find('-')

  desiredvalues = np.array(dimt, dtype=bool)
  desiredvaluesi = np.array(dimt, dtype=bool)
  desiredvaluese = np.array(dimt, dtype=bool)

  if not isper == -1:
# Months list given as a period [MMi]-[MMf]
    print '      submns: There is a period of months "' + values + '"'
    imn = int(values.split('-')[0])
    emn = int(values.split('-')[1])
    desiredvaluesi = np.array(realdates[:,1] >= imn)
    desiredvaluese = np.array(realdates[:,1] <= emn)

    desiredvalues = desiredvaluesi*desiredvaluese

  else:
    mns = values.split(':')
    nmn = 1
    for imn in range(len(mns)):
      desiredvaluesi = np.array(realdates[:,1] == int(mns[imn]))
      if nmn == 1:
        desiredvalues = desiredvaluesi
      else:
        desiredvalues = desiredvaluesi+desiredvalues
      nmn = nmn + 1

  Ndesiredvalues = len(realdates[desiredvalues])

  print '    submns: N values: ', Ndesiredvalues
  if Ndesiredvalues == 0:
    print errormsg
    print '    submns: No values found for "' + values + '"!!'
    ncf.close()
    quit(-1)

# Variable values (assuming time as first dimension)
##
  var = ncf.variables[varn]
#  varvals = var[:]
  vardims = var.shape
  varshape = len(vardims)
#  print '      submns: Shape of data: ',varshape, ':' , vardims

  if varshape == 1:
    vardesiredvalues = np.arange(Ndesiredvalues)
    vardesiredvalues = var[desiredvalues]
  elif varshape == 2:
    vardesiredvalues = np.arange(Ndesiredvalues*vardims[1]).reshape(Ndesiredvalues,vardims[1])
    vardesiredvalues = var[desiredvalues, :]
  elif varshape == 3:
    vardesiredvalues = np.arange(Ndesiredvalues*vardims[1]*vardims[2]).reshape(Ndesiredvalues,                        \
      vardims[1],vardims[2])
    vardesiredvalues = var[desiredvalues, :, :]
  elif varshape == 4:
    vardesiredvalues = np.arange(Ndesiredvalues*vardims[1]*vardims[2]*vardims[3]).reshape(Ndesiredvalues,             \
      vardims[1],vardims[2],vardims[3])
    vardesiredvalues = var[desiredvalues, :, :, :]
  elif varshape == 5:
    vardesiredvalues = np.arange(Ndesiredvalues*vardims[1]*vardims[2]*vardims[3]*vardims[4]).reshape(Ndesiredvalues,  \
      vardims[1],vardims[2],vardims[3],vardims[4])
    vardesiredvalues = var[desiredvalues, :, :, :, :]
  elif varshape == 6:
    vardesiredvalues = np.arange(Ndesiredvalues*vardims[1]*vardims[2]*vardims[3]*vardims[4]*                          \
      vardims[5]).reshape(Ndesiredvalues,vardims[1],vardims[2],vardims[3],vardims[4],vardims[5])
    vardesiredvalues = var[desiredvalues, :, :, :, :, :]
  else:
     print errormsg
     print '    submns: ', varshape, ' shape of matrix not prepared !'
     ncf.close()
     quit(-1)

#  print '  shape of desired values: ', vardesiredvalues.shape
# Creation of file
##

  ncfo = NetCDFFile( ofile, 'w')

  vardims = var.dimensions
  vartype = var.dtype
  varattr = var.ncattrs()

# Checking dimensions
##
  newdims = ncfo.dimensions
  for rdim in vardims:
      if not searchInlist(newdims, rdim):
          if not rdim == 'time':
              print '      submns: Adding dimension ' + rdim
              ncfo.sync()
              ncf.close()
              ncfo.close()
              fdimadd(ncfile + ':' + rdim, ofile)
              ncf = NetCDFFile(ncfile,'r')
              ncfo = NetCDFFile(ofile,'a')
          else:
              ncfo.createDimension('time', None)

# Checking fill value
## 
  if searchInlist(varattr, '_FillValue'):
      varfil = var._FillValue
  else:
      varfil = False

  newvar = ncfo.createVariable(varn, vartype, vardims, fill_value=varfil)

  if not varshape == 0:
    newvar[:] = vardesiredvalues

  newvar = ncfo.variables[varn]
  for attr in varattr:
      newvarattrs = newvar.ncattrs()
      attrv = var.getncattr(attr)
      if not searchInlist(newvarattrs, attr):     
          newvar.setncattr(attr, attrv)

  vardims = times.dimensions
  vartype = times.dtype
  varattr = times.ncattrs()

  newvar = ncfo.createVariable('time', vartype, vardims, fill_value=varfil)
  newvar = ncfo.variables['time']
  newvar[:] = timevals[desiredvalues]

  ncf.close()
  ncfo.sync()
  ncfo.close()
  fattradd('time', ncfile + ':time', ofile)
  fvaradd(ncfile + ':lon', ofile)
  fvaradd(ncfile + ':lat', ofile)

  ncfo = NetCDFFile(ofile,'a')
  newvar = ncfo.variables['time']
  newvarattrs = newvar.ncattrs()
  if searchInlist(newvarattrs, 'bounds'):
      if newvar.getncattr('bounds') == 'time_bnds':
          ncf = NetCDFFile(ncfile,'r')
          tbnds = ncf.variables['time_bnds']
          vardims = tbnds.dimensions
          vartype = tbnds.dtype
          varattr = tbnds.ncattrs()
          ncfo.createDimension('bnds', 2)
          newvar = ncfo.createVariable('time_bnds', vartype, vardims, fill_value=varfil)
          newvar[:] = tbnds[desiredvalues,:]

          ncf.close()
          ncfo.sync()
          ncfo.close()
          fattradd('time_bnds', ncfile + ':time_bnds', ofile)
      else:
          ncfo.close()
  else:
      ncfo.close()

  fgaddattr(ncfile, ofile)

  print '      submns: File "' + ofile + '" with a subset of ' + values + ' has been created'

class statsValWeigthed(object):
  """Weigthed Statistics class providing:
  vals = values (can be a matrix)
  wgs = weights (can be a matrix)
  self.meanv: mean weigthed value
  self.mean2v: mean quadratic weigthed value
  self.stdv: weigthed standard deviation
  self.Nokvalue non None values of a list of values 
  self.meanwgt: mean of the weigths
  self.mean2wgt: cuadratic mean of the weigths
  self.stdwgt: standard deviation of the weigths
  """

  def __init__(self, vals, wgs):
    import math

    if vals is None:
      self.Nv = None
      self.meanv = None
      self.mean2v = None
      self.stdv = None
      self.Nokvalues = None
      self.meanwgt = None
      self.mean2wgt = None
      self.stdwgt = None
    else:
      values = vals.flat 
      weights = wgs.flat
      self.Nv=len(values)
      self.meanv=0.
      self.mean2v=0.
      self.stdv=0.
      self.meanwgt = 0.
      self.mean2wgt = 0.
      self.stdwgt = 0.
      self.Nokvalues = 0

      for inum in range(self.Nv):
        if not values[inum] is None:
          self.Nokvalues = self.Nokvalues + 1
          self.meanv = self.meanv+values[inum]*weights[inum]
          self.mean2v = self.mean2v+values[inum]*weights[inum]*values[inum]
          self.meanwgt = self.meanwgt+weights[inum]
          self.mean2wgt = self.mean2wgt+weights[inum]*weights[inum]

      self.meanv = self.meanv/float(self.meanwgt)
      self.mean2v = self.mean2v/float(self.meanwgt)
      self.stdv = math.sqrt(self.mean2v-self.meanv*self.meanv)
      self.meanwgt = self.meanwgt/float(self.Nokvalues)
      self.mean2wgt = self.mean2wgt/float(self.Nokvalues)
      self.stdwgt = math.sqrt(self.mean2wgt-self.meanwgt*self.meanwgt)

class stats2Val(object):
  """two variables Statistics class providing:
  vals1 = variable 1
  vals2 = variable 2
  power = power of the polynomial fitting to apply between both variables
  self.min[var], self.max[var], self.mean[var], self.mean2[var], self.std[var] of 
    [var] = var1+var2[v1Av2], var1-var2[v1Sv2], var1/var2[v1Dv2], var1*var2[v1Pv2]
  self.Nokvalues1: number of correct values of variable 1
  self.Nokvalues2: number of correct values of variable 2
  self.Nokvalues12: number of correct coincident values of variable 1 and variable 2
  self.mae=mean(abs(var1-var2)) 
  self.rmse=sqrt((var1-var2)**2) 
  self.correlation (and p-value) 
  self.linRegress: linear regression [trend, intercept, regression coefficient, p_value, standard error]
  self.polRegress: polinomial Regresion  of degree [power] [coef**[power], coef**[power-1], ...., coef**0]
  """

  def __init__(self, vals1, vals2, power):
    import math
    import numpy as np
    from scipy import stats as sts

    if vals1 is None:
      self.Nv = None
      self.Nokvalues1 = None
      self.Nokvalues2 = None
      self.Nokvalues12 = None
      self.NDvalNone = None
      self.minv1Av2 = None
      self.maxv1Av2 = None
      self.meanv1Av2 = None
      self.mean2v1Av2 = None
      self.stdv1Av2 = None
      self.minv1Sv2 = None
      self.maxv1Sv2 = None
      self.meanv1Sv2 = None
      self.mean2v1Sv2 = None
      self.stdv1Sv2 = None
      self.minv1Dv2 = None
      self.maxv1Dv2 = None
      self.meanv1Dv2 = None
      self.mean2v1Dv2 = None
      self.stdv1Dv2 = None
      self.minv1Pv2 = None
      self.maxv1Pv2 = None
      self.meanv1Pv2 = None
      self.mean2v1Pv2 = None
      self.stdv1Pv2 = None
      self.mae = None
      self.rmse = None
      self.corr = None
      self.linRegress = None
      self.polRegress = None
      self.polRegressResidual = None
      self.polRegressRes = None
      self.polRegressSingVal = None
    else:
      values1 = vals1.flat 
      values2 = vals2.flat 

      if not len(values1) == len(values2):
        print errormsg
        print '    stats2Val: lengths of variables differ!! Lvar1: ', len(values1), ' Lvar2: ',len(values2),' statistics between them can not be computed!'
        quit(-1)

      self.Nv=len(values1)
      self.minv1Av2=10000000000.
      self.maxv1Av2=-self.minv1Av2
      self.meanv1Av2=0.
      self.mean2v1Av2=0.
      self.stdv1Av2=0.
      self.minv1Sv2=self.minv1Av2
      self.maxv1Sv2=-self.minv1Av2
      self.meanv1Sv2=0.
      self.mean2v1Sv2=0.
      self.stdv1Sv2=0.
      self.minv1Dv2=self.minv1Av2
      self.maxv1Dv2=-self.minv1Av2
      self.meanv1Dv2=0.
      self.mean2v1Dv2=0.
      self.stdv1Dv2=0.
      self.minv1Pv2=self.minv1Av2
      self.maxv1Pv2=-self.minv1Av2
      self.meanv1Pv2=0.
      self.mean2v1Pv2=0.
      self.stdv1Pv2=0.
      self.mae = 0.
      self.rmse = 0.
      self.corr = np.array([0., 0.])
      self.linRegress = np.zeros(5, float)
      self.polRegress = np.zeros(power+1, float)
      self.polRegressResidual = 0.
      self.polRegressSingVal = np.zeros(power+1, float)

# v1 [+ / - / / / *] v2
##
      self.Nokvalues1 = 0
      self.Nokvalues2 = 0
      self.Nokvalues12 = 0
      self.NDvalNone = 0
      for inum in range(self.Nv):
        if not values1[inum] is None:
          self.Nokvalues1 = self.Nokvalues1 + 1
        if not values2[inum] is None:
          self.Nokvalues2 = self.Nokvalues2 + 1
        if not values1[inum] is None and not values2[inum] is None:
          self.Nokvalues12 = self.Nokvalues12 + 1
          Aval = values1[inum] + values2[inum]
          Sval = values1[inum] - values2[inum]
          Pval = values1[inum] * values2[inum]
          if math.isinf(values1[inum] / values2[inum]) or math.isnan(values1[inum] / values2[inum]):
            if self.NDvalNone < 1:
               print warnmsg
               print '      stats2Val: val1/val2 inf or Nan!!!!'
            Dval = None
            self.NDvalNone = self.NDvalNone + 1
          else:
            Dval = values1[inum] / values2[inum]

          self.mae = self.mae + abs(Sval)
          self.rmse = self.rmse + Sval**2

          if Aval < self.minv1Av2:
            self.minv1Av2 = Aval
          if Aval > self.maxv1Av2:
            self.maxv1Av2 = Aval
          if Sval < self.minv1Sv2:
            self.minv1Sv2 = Sval
          if Sval > self.maxv1Sv2:
            self.maxv1Sv2 = Sval
          if not Dval is None and Dval < self.minv1Dv2:
            self.minv1Dv2 = Dval
          if not Dval is None and  Dval > self.maxv1Dv2:
            self.maxv1Dv2 = Dval
          if Pval < self.minv1Pv2:
            self.minv1Pv2 = Pval
          if Pval > self.maxv1Pv2:
            self.maxv1Pv2 = Pval

          self.meanv1Av2 = self.meanv1Av2+Aval
          self.mean2v1Av2 = self.mean2v1Av2+Aval*Aval
          self.meanv1Sv2 = self.meanv1Sv2+Sval
          self.mean2v1Sv2 = self.mean2v1Sv2+Sval*Sval
          if not Dval is None:
            self.meanv1Dv2 = self.meanv1Dv2+Dval
            self.mean2v1Dv2 = self.mean2v1Dv2+Dval*Dval
          self.meanv1Pv2 = self.meanv1Pv2+Pval
          self.mean2v1Pv2 = self.mean2v1Pv2+Pval*Pval

##      print 'Nokvalues1: ', self.Nokvalues1, 'Nokvalues2: ', self.Nokvalues2, 'Nokvalues12: ', float(self.Nokvalues12), 'NDvalNone: ',self.NDvalNone
      self.meanv1Av2 = self.meanv1Av2/float(self.Nokvalues12)
      self.mean2v1Av2 = self.mean2v1Av2/float(self.Nokvalues12)
      self.stdv1Av2 = math.sqrt(self.mean2v1Av2-self.meanv1Av2*self.meanv1Av2)
      self.meanv1Sv2 = self.meanv1Sv2/float(self.Nokvalues12)
      self.mean2v1Sv2 = self.mean2v1Sv2/float(self.Nokvalues12)
      self.stdv1Sv2 = math.sqrt(self.mean2v1Sv2-self.meanv1Sv2*self.meanv1Sv2)
      if self.Nokvalues12 - self.NDvalNone == 0:
          self.meanv1Dv2 = None
          self.mean2v1Dv2 = None
          self.stdv1Dv2 = None
          print warnmsg
          print '      stats2Val: all values of val1/val2 are None!'
      else:
          self.meanv1Dv2 = self.meanv1Dv2/(float(self.Nokvalues12 - self.NDvalNone))
          self.mean2v1Dv2 = self.mean2v1Dv2/(float(self.Nokvalues12 - self.NDvalNone))
          self.stdv1Dv2 = math.sqrt(self.mean2v1Dv2-self.meanv1Dv2*self.meanv1Dv2)
      self.meanv1Pv2 = self.meanv1Pv2/float(self.Nokvalues12)
      self.mean2v1Pv2 = self.mean2v1Pv2/float(self.Nokvalues12)
      self.stdv1Pv2 = math.sqrt(self.mean2v1Pv2-self.meanv1Pv2*self.meanv1Pv2)

      self.mae = self.mae/self.Nokvalues12
      self.rmse = math.sqrt(self.rmse/self.Nokvalues12)

      self.corr = sts.pearsonr(values1, values2)

#      linresMat = np.array([ vals1, np.ones(len(vals1)) ]).T
      self.linRegress[0], self.linRegress[1], self.linRegress[2], self.linRegress[3], self.linRegress[4] = sts.linregress(vals1, vals2)
#      self.linRegress[0], self.linRegress[1], self.linRegress[2], self.linRegress[3] = sts.linregress(linresMat, vals2)

      polyfitvals=np.polyfit(values1, values2, power, full = True)

      self.polRegress = polyfitvals[0]
      self.polRegressRes = polyfitvals[1]
      self.polRegressSingVal = polyfitvals[3]

class statsVal(object):
  """Statistics class providing
  vals = variable
  self.Nv = number of values
  self.minv = minimum value
  self.maxv = maximum value
  self.meanv = mean value
  self.mean2v = cuadratic mean value
  self.stdv = standard deviation value
  self.Nokvalues = number of correct values of variable
  self.quantilesv = quantiles (%5 bins) of the variable
  """

  def __init__(self, vals):
    import math

    if vals is None:
      self.Nv = None
      self.minv = None
      self.maxv = None
      self.meanv = None
      self.mean2v = None
      self.stdv = None
      self.Nokvalues = None
      self.quantilesv = None
    else:
      values = vals.flat 
      self.Nv=len(values)
      self.minv=10000000000.
      self.maxv=-self.minv
      self.meanv=0.
      self.mean2v=0.
      self.stdv=0.

      sortedvalues = sorted(values)

      self.Nokvalues = 0
      for inum in range(self.Nv):
        if not values[inum] is None:
          self.Nokvalues = self.Nokvalues + 1
          if values[inum] < self.minv:
            self.minv = values[inum]
          if values[inum] > self.maxv:
            self.maxv = values[inum]

          self.meanv = self.meanv+values[inum]
          self.mean2v = self.mean2v+values[inum]*values[inum]

      self.meanv = self.meanv/float(self.Nokvalues)
      self.mean2v = self.mean2v/float(self.Nokvalues)
      self.stdv = math.sqrt(self.mean2v-self.meanv*self.meanv)
      self.quantilesv = []
      for iq in range(20):
        self.quantilesv.append(sortedvalues[int((self.Nv-1)*iq/20)])

      self.quantilesv.append(sortedvalues[self.Nv-1])
      self.medianv = self.quantilesv[10]

def spacemean(ncfile, varn):
  """ Function to retrieve a space mean series from a multidimensional variable of a file
  ncfile = netCDF file name
  varn = variable name
  """
  import datetime as dt
  import calendar as cal

  ofile = 'spacemean_' + varn + '.nc'
  varfil=1.e20
  statsn = ['min', 'max', 'mean', 'mean2', 'stdv', 'meanwgt', 'mean2wgt', 'stdvwgt', 'quant']
  statslongn = ['minimum', 'maximum', 'mean', 'quadratic mean', 'standard deviation', 'weigthed mean',  'weigthed quadratic mean', 'weigthed standard deviation', 'quantiles']

  ncf = NetCDFFile(ncfile,'r')

  if not ncf.variables.has_key(varn):
    print errormsg
    print '    spacemean: File "' + ncfile + '" does not have variable "' + varn + '" !!!!'
    print errormsg
    ncf.close()
    quit(-1)

  times = ncf.variables['time']
  timevals = times[:]

  attvar = times.ncattrs()
  if not searchInlist(attvar, 'units'):
    print errormsg
    print '    spacemean: "time" does not have attribute: "units"'
    ncf.close()
    quit(-1)
  else:
    units = times.getncattr('units')
  
  txtunits = units.split(' ')
  tunits = txtunits[0]
  Srefdate = txtunits[len(txtunits) - 1]
# Does reference date contain a time value [YYYY]-[MM]-[DD] [HH]:[MI]:[SS]
##
  timeval = Srefdate.find(':')

  if not timeval == -1:
    print '      spacemean: refdate with time!'
    refdate = datetimeStr_datetime(txtunits[len(txtunits) - 2] + '_' + Srefdate)
  else:
    refdate = dateStr_date(Srefdate)

  dimt = len(timevals)
  realdates = np.zeros((dimt, 6), dtype=int)
  print realdates.shape

## Not in timedelta
#  if tunits == 'years':
#    for it in range(dimt):
#      realdate = refdate + dt.timedelta(years=float(times[it]))
#      realdates[it] = int(realdate.year)
#  elif tunits == 'months':
#    for it in range(dimt):
#      realdate = refdate + dt.timedelta(months=float(times[it]))
#      realdates[it] = int(realdate.year)
  if tunits == 'weeks':
    for it in range(dimt):
      realdate = refdate + dt.timedelta(weeks=float(times[it]))
      realdates[it,0] = int(realdate.year)
      realdates[it,1] = int(realdate.month)
      realdates[it,2] = int(realdate.day)
      realdates[it,3] = int(realdate.hour)
      realdates[it,4] = int(realdate.second)
      realdates[it,5] = int(realdate.minute)
  elif tunits == 'days':
    for it in range(dimt):
      realdate = refdate + dt.timedelta(days=float(times[it]))
      realdates[it,0] = int(realdate.year)
      realdates[it,1] = int(realdate.month)
      realdates[it,2] = int(realdate.day)
      realdates[it,3] = int(realdate.hour)
      realdates[it,4] = int(realdate.second)
      realdates[it,5] = int(realdate.minute)
  elif tunits == 'hours':
    for it in range(dimt):
      realdate = refdate + dt.timedelta(hours=float(times[it]))
      realdates[it,0] = int(realdate.year)
      realdates[it,1] = int(realdate.month)
      realdates[it,2] = int(realdate.day)
      realdates[it,3] = int(realdate.hour)
      realdates[it,4] = int(realdate.second)
      realdates[it,5] = int(realdate.minute)
  elif tunits == 'minutes':
    for it in range(dimt):
      realdate = refdate + dt.timedelta(minutes=float(times[it]))
      realdates[it,0] = int(realdate.year)
      realdates[it,1] = int(realdate.month)
      realdates[it,2] = int(realdate.day)
      realdates[it,3] = int(realdate.hour)
      realdates[it,4] = int(realdate.second)
      realdates[it,5] = int(realdate.minute)
  elif tunits == 'seconds':
    for it in range(dimt):
      realdate = refdate + dt.timedelta(seconds=float(times[it]))
      realdates[it,0] = int(realdate.year)
      realdates[it,1] = int(realdate.month)
      realdates[it,2] = int(realdate.day)
      realdates[it,3] = int(realdate.hour)
      realdates[it,4] = int(realdate.second)
      realdates[it,5] = int(realdate.minute)
  elif tunits == 'milliseconds':
    for it in range(dimt):
      realdate = refdate + dt.timedelta(milliseconds=float(times[it]))
      realdates[it,0] = int(realdate.year)
      realdates[it,1] = int(realdate.month)
      realdates[it,2] = int(realdate.day)
      realdates[it,3] = int(realdate.hour)
      realdates[it,4] = int(realdate.second)
      realdates[it,5] = int(realdate.minute)
  else:
    print errormsg
    print '    spacemean: time units "' + tunits + '" not ready!!!!'
    ncf.close()
    quit(-1)

# Variable values (assuming time as first dimension)
##
  var = ncf.variables[varn]
  varinf = variable_inf(var)
  if searchInlist(varinf.attributes, 'scale_factor'):
    scalefact = var.getncattr('scale_factor')
    print '      spacemean: data has a scale factor of ', scalefact
  else:
    scalefact = 1.

  if searchInlist(varinf.attributes, 'add_offset'):
    offset = var.getncattr('add_offset')
    print '      spacemean: data has an offset of ', offset
  else:
    offset = 0.

#  varvals = var[:]
  vardims = var.shape
  varshape = len(vardims)
  vardimns = var.dimensions
##  print '      spacemean: Shape of data: ',varshape, ':' , vardims
  dimt=vardims[0]
  vardnames = []

# Spatial average spatial weigthed
## 
  lonvar = ncf.variables['lon']
  latvar = ncf.variables['lat']

  if not len(lonvar.shape) == 2:
    lonv = lonvar[:]
    latv = latvar[:]
    dx=lonvar.shape[0]
    dy=latvar.shape[0]
    lonval = np.zeros((dy, dx), dtype=float)    
    latval = np.zeros((dy, dx), dtype=float)    
    for iy in range(dy):
      lonval[iy,:] = lonv
    for ix in range(dx):
      latval[:,ix] = latv
  else:
    lonval = lonvar[:]
    latval = latvar[:]

  weightsv = abs(np.cos(latval*np.pi/180.))

  if varinf.Ndims == 1:
    print errormsg
    print '    spacemean: You can not compute a space mean for a ', varinf.Ndims, 'D var!!!'
    ncf.close()
    quit(-1)
  elif varinf.Ndims== 2:
    print errormsg
    print '    spacemean: You can not compute a space mean for a ', varinf.Ndims, 'D var!!!'
    ncf.close()
    quit(-1)
  elif varinf.Ndims == 3:
    varstats = np.ones((dimt, 8), dtype=float)
    varstats = varstats*varfil
    varquant = np.ones((dimt, 21), dtype=float)
    varquant = varquant*varfil
    vardnames.append(vardimns[0])
    dy=vardims[1]
    dx=vardims[2]
    for it in range(dimt):
      if not varinf.FillValue is None:
        varvl = var[it,:,:]
        varval = np.where(varvl == varinf.FillValue, None, varvl)
      else:
        varval = var[it,:,:]

      varst = statsVal(varval)
      varstwgt  = statsValWeigthed(varval, weightsv)
      varstats[it,0] = varst.minv
      varstats[it,1] = varst.maxv
      varstats[it,2] = varst.meanv
      varstats[it,3] = varst.mean2v
      varstats[it,4] = varst.stdv
      varquant[it,:] = varst.quantilesv
      varstats[it,5] = varstwgt.meanv
      varstats[it,6] = varstwgt.mean2v
      varstats[it,7] = varstwgt.stdv
  elif varshape == 4:
    varstats = np.ones((dimt, vardims[1], 8), dtype=float)
    varstats = varstats*varfil
    varquant = np.ones((dimt, vardims[1], 21), dtype=float)
    varquant = varquant*varfil
    vardimnames = (str(vardimns[0]), str(vardimns[1]))
    vardnames.append(vardimns[0])
    vardnames.append(vardimns[1])
    dy=vardims[2]
    dx=vardims[3]
    for it in range(dimt):
      for ik in range(vardims[1]):
        if not varinf.FillValue is None:
          varvl = var[it,ik,:,:]
          varval = np.where(varvl == varinf.FillValue, None, varvl)
        else:
          varval = var[it,ik,:,:]

        varst = statsVal(varval)
        varstwgt  = statsValWeigthed(varval, weightsv)
        varstats[it,ik,0] = varst.minv
        varstats[it,ik,1] = varst.maxv
        varstats[it,ik,2] = varst.meanv
        varstats[it,ik,3] = varst.mean2v
        varstats[it,ik,4] = varst.stdv
        varquant[it,ik,:] = varst.quantilesv
        varstats[it,ik,5] = varstwgt.meanv
        varstats[it,ik,6] = varstwgt.mean2v
        varstats[it,ik,7] = varstwgt.stdv
  elif varshape == 5:
    varstats = np.ones((dimt, vardims[1], vardims[2], 8), dtype=float)
    varstats = varstats*varfil
    varquant = np.ones((dimt,  vardims[1], vardims[2], 21), dtype=float)
    varquant = varquant*varfil
    vardnames.append(vardimns[0])
    vardnames.append(vardimns[1])
    vardnames.append(vardimns[2])
    dy=vardims[3]
    dx=vardims[4]
    for it in range(dimt):
      for ik in range(vardims[1]):
        for il in range(vardims[2]):
          if not varinf.FillValue is None:
            varvl = var[it,ik,il,:,:]
            varval = np.where(varvl == varinf.FillValue, None, varvl)
          else:
            varval = var[it,ik,il,:,:]

          varst = statsVal(varval)
          varstwgt  = statsValWeigthed(varval, weightsv)
          varstats[it,ik,il,0] = varst.minv
          varstats[it,ik,il,1] = varst.maxv
          varstats[it,ik,il,2] = varst.meanv
          varstats[it,ik,il,3] = varst.mean2v
          varstats[it,ik,il,4] = varst.stdv
          varquant[it,ik,il,:] = varst.quantilesv
          varstats[it,ik,il,5] = varstwgt.meanv
          varstats[it,ik,il,6] = varstwgt.mean2v
          varstats[it,ik,il,7] = varstwgt.stdv
  else:
     print errormsg
     print '    spacemean: ', varshape, ' shape of matrix not prepared !'
     ncf.close()
     quit(-1)

  vardimnames = tuple(vardnames)
#  print '  shape of desired values: ', vardesiredvalues.shape
# Creation of file
##

  ncfo = NetCDFFile( ofile, 'w')

  vartype = var.dtype
  varattr = var.ncattrs()

# Checking dimensions
##
  newdims = ncfo.dimensions
  Nvardnames = len(vardimnames)
  for idim in range(Nvardnames):
      rdim = vardimnames[idim] 
      if not searchInlist(newdims, rdim):
          if not rdim == 'time':
              print '      spacemean: Adding dimension ' + rdim
              ncfo.sync()
              ncf.close()
              ncfo.close()
              fdimadd(ncfile + ':' + rdim, ofile)
              ncf = NetCDFFile(ncfile,'r')
              ncfo = NetCDFFile(ofile,'a')
          else:
              ncfo.createDimension('time', None)

# Checking fill value
## 
  if searchInlist(varattr, '_FillValue'):
      varfil = var._FillValue
  else:
      varfil = False

  Nstats = len(statsn)
  for ist in range(Nstats):
    if statsn[ist] == 'quant':
      print ist, statsn[ist]##, ': ', varquant[int(dimt/2),10]
      newdim = ncfo.createDimension('quant', 21)
      newvardnames = list(vardnames)
      newvardnames.append('quant')
      newvar = ncfo.createVariable(varn + statsn[ist], vartype, tuple(newvardnames), fill_value=varfil)
      newvar[:] = varquant
    else:
      print ist, statsn[ist]##, ': ', varstats[int(dimt/2),ist]
      newvar = ncfo.createVariable(varn + statsn[ist], vartype, vardimnames, fill_value=varfil)
      newvar[:] = varstats[:,ist]

    newvar = ncfo.variables[varn + statsn[ist]]
    for attr in varattr:
         newvarattrs = newvar.ncattrs()
         attrv = var.getncattr(attr)
         if not searchInlist(newvarattrs, attr):
             if attr == 'scale_factor' or attr == 'add_offset':
                 if not statsn[ist] == 'stdv' and not statsn[ist] == 'stdvwgt':
                     newvar.setncattr(attr, attrv)
             else:
                 newvar.setncattr(attr, attrv)
    
    newvar.setncattr('cell_methods', 'space ' + statslongn[ist] + ' all domain ' + str(dy) + 'x' + str(dx))

##  print '  Adding time variable'
  vartdims = times.dimensions
  vartype = times.dtype
  varattr = times.ncattrs()

  newvar = ncfo.createVariable('time', vartype, vartdims, fill_value=varfil)
  newvar = ncfo.variables['time']
  newvar[:] = times

  ncf.close()
  ncfo.sync()
  ncfo.close()
  fattradd('time', ncfile + ':time', ofile)
  fvaradd(ncfile + ':lon', ofile)
  fvaradd(ncfile + ':lat', ofile)

  ncfo = NetCDFFile(ofile,'a')
  newvar = ncfo.variables['time']
  newvarattrs = newvar.ncattrs()
  if searchInlist(newvarattrs, 'bounds'):
      if newvar.getncattr('bounds') == 'time_bnds':
          ncf = NetCDFFile(ncfile,'r')
          tbnds = ncf.variables['time_bnds']
          vardims = tbnds.dimensions
          vartype = tbnds.dtype
          varattr = tbnds.ncattrs()
          ncfo.createDimension('bnds', 2)
          newvar = ncfo.createVariable('time_bnds', vartype, vardims, fill_value=varfil)
          newvar[:] = tbnds[:]

          ncf.close()
          ncfo.sync()
          ncfo.close()
          fattradd('time_bnds', ncfile + ':time_bnds', ofile)
      else:
          ncfo.close()
  else:
      ncfo.close()

  fgaddattr(ncfile, ofile)

  print '    spacemean: File "' + ofile + '" as space mean of "' + varn + '" has been created'

def timemean(values, ncfile, varn):
  """ Function to retrieve a time mean series from a multidimensional variable of a file
  values = power of the polynomial fitting with time to be applied
  ncfile = netCDF file name
  varn = variable name
  """
  import datetime as dt
  import calendar as cal

  powerv=int(values)
  ofile = 'timemean_' + varn + '.nc'
  varfil=1.e20
  statsn = ['min', 'max', 'mean', 'mean2', 'stdv', 'quant','linregress','polregress']
  statslongn = ['minimum', 'maximum', 'mean', 'quadratic mean', 'standard deviation', 'quantiles', \
    'linear regression', 'polynomial regression']

  ncf = NetCDFFile(ncfile,'r')

  if not ncf.variables.has_key(varn):
    print errormsg
    print '    timemean: File  "' + ncfile + '" does not have variable ' + varn + ' !!!!'
    print errormsg
    ncf.close()
    quit(-1)

  times = ncf.variables['time']
  timevals = times[:]

  attvar = times.ncattrs()
  if not searchInlist(attvar, 'units'):
    print errormsg
    print '    timemean: "time" does not have attribute: "units"'
    ncf.close()
    quit(-1)
  else:
    units = times.getncattr('units')
  
  txtunits = units.split(' ')
  tunits = txtunits[0]
  Srefdate = txtunits[len(txtunits) - 1]
# Does reference date contain a time value [YYYY]-[MM]-[DD] [HH]:[MI]:[SS]
##
  timeval = Srefdate.find(':')

  if not timeval == -1:
    print '      timemean: refdate with time!'
    refdate = datetimeStr_datetime(txtunits[len(txtunits) - 2] + '_' + Srefdate)
  else:
    refdate = dateStr_date(Srefdate)

  dimt = len(timevals)
  realdates = np.zeros((dimt, 6), dtype=int)
  print realdates.shape

## Not in timedelta
#  if tunits == 'years':
#    for it in range(dimt):
#      realdate = refdate + dt.timedelta(years=float(times[it]))
#      realdates[it] = int(realdate.year)
#  elif tunits == 'months':
#    for it in range(dimt):
#      realdate = refdate + dt.timedelta(months=float(times[it]))
#      realdates[it] = int(realdate.year)
  if tunits == 'weeks':
    for it in range(dimt):
      realdate = refdate + dt.timedelta(weeks=float(times[it]))
      realdates[it,0] = int(realdate.year)
      realdates[it,1] = int(realdate.month)
      realdates[it,2] = int(realdate.day)
      realdates[it,3] = int(realdate.hour)
      realdates[it,4] = int(realdate.second)
      realdates[it,5] = int(realdate.minute)
  elif tunits == 'days':
    for it in range(dimt):
      realdate = refdate + dt.timedelta(days=float(times[it]))
      realdates[it,0] = int(realdate.year)
      realdates[it,1] = int(realdate.month)
      realdates[it,2] = int(realdate.day)
      realdates[it,3] = int(realdate.hour)
      realdates[it,4] = int(realdate.second)
      realdates[it,5] = int(realdate.minute)
  elif tunits == 'hours':
    for it in range(dimt):
      realdate = refdate + dt.timedelta(hours=float(times[it]))
      realdates[it,0] = int(realdate.year)
      realdates[it,1] = int(realdate.month)
      realdates[it,2] = int(realdate.day)
      realdates[it,3] = int(realdate.hour)
      realdates[it,4] = int(realdate.second)
      realdates[it,5] = int(realdate.minute)
  elif tunits == 'minutes':
    for it in range(dimt):
      realdate = refdate + dt.timedelta(minutes=float(times[it]))
      realdates[it,0] = int(realdate.year)
      realdates[it,1] = int(realdate.month)
      realdates[it,2] = int(realdate.day)
      realdates[it,3] = int(realdate.hour)
      realdates[it,4] = int(realdate.second)
      realdates[it,5] = int(realdate.minute)
  elif tunits == 'seconds':
    for it in range(dimt):
      realdate = refdate + dt.timedelta(seconds=float(times[it]))
      realdates[it,0] = int(realdate.year)
      realdates[it,1] = int(realdate.month)
      realdates[it,2] = int(realdate.day)
      realdates[it,3] = int(realdate.hour)
      realdates[it,4] = int(realdate.second)
      realdates[it,5] = int(realdate.minute)
  elif tunits == 'milliseconds':
    for it in range(dimt):
      realdate = refdate + dt.timedelta(milliseconds=float(times[it]))
      realdates[it,0] = int(realdate.year)
      realdates[it,1] = int(realdate.month)
      realdates[it,2] = int(realdate.day)
      realdates[it,3] = int(realdate.hour)
      realdates[it,4] = int(realdate.second)
      realdates[it,5] = int(realdate.minute)
  else:
    print errormsg
    print '    timemean: time units "' + tunits + '" not ready!!!!'
    ncf.close()
    quit(-1)

  timesv = (realdates[:,0] - realdates[0,0])*12 + realdates[:,1]

# Variable values (assuming time as first dimension)
##
  var = ncf.variables[varn]
  varinf = variable_inf(var)

#  varvals = var[:]
  vardims = varinf.dims
  varshape = varinf.Ndims
  vardimns = var.dimensions

##  print '     timemean: Shape of data: ',varshape, ':' , vardims
  dimt=vardims[0]
  vardnames = []

  if varshape == 1:
    print errormsg
    print '    timemean: You can not compute a time mean for a ', varshape, 'D var!!!'
    ncf.close()    
    quit(-1)
  elif varshape == 2:
    dx=vardims[1]
    varstats = np.ones((dx, 5), dtype=float)
    varstats = varstats*varfil
    varquant = np.ones((dx, 21), dtype=float)
    varquant = varquant*varfil
    varlregress = np.ones((dx, 5), dtype=float)
    varlregress = varlregress*varfil
    varpregress = np.ones((dx, powerv+1), dtype=float)
    varpregress = varpregress*varfil
    varpregressres = np.ones((dx), dtype=float)
    varpregressres = varpregressres*varfil
    varpregresssingval = varpregress
    varpregresssingval = varpregresssingval*varfil
    vardnames.append(vardimns[1])
    for ix in range(dx):
      if not varinf.FillValue is None:
        varvl = var[:,ix]
        varval = np.where(varvl == varinf.FillValue, None, varvl)
      else:
        varval = var[:,ix]
      
      varst = statsVal(varval)
      var2st = stats2Val(timesv, varval, powerv)
      varstats[ix,0] = varst.minv
      varstats[ix,1] = varst.maxv
      varstats[ix,2] = varst.meanv
      varstats[ix,3] = varst.mean2v
      varstats[ix,4] = varst.stdv
      varquant[ix,:] = varst.quantilesv
      varlregress[ix,:] = var2st.linRegress
      varpregress[ix,:] = var2st.polRegress
      varpregressres[ix] = var2st.polRegressRes
      varpregresssingval[ix,:] = var2st.polRegressSingVal
  elif varshape == 3:
    dy=vardims[1]
    dx=vardims[2]
    varstats = np.ones((dy, dx, 5), dtype=float)
    varstats = varstats*varfil
    varquant = np.ones((dy, dx, 21), dtype=float)
    varquant = varquant*varfil
    varlregress = np.ones((dy, dx, 5), dtype=float)
    varlregress = varlregress*varfil
    varpregress = np.ones((dy, dx, powerv+1), dtype=float)
    varpregress = varpregress*varfil
    varpregressres = np.ones((dy,dx), dtype=float)
    varpregressres = varpregressres*varfil
    varpregresssingval = varpregress
    varpregresssingval = varpregresssingval*varfil
    vardnames.append(vardimns[1])
    vardnames.append(vardimns[2])
    for iy in range(dy):
      for ix in range(dx):
        if not varinf.FillValue is None:
          varvl = var[:,iy,ix]
          varval = np.where(varvl == varinf.FillValue, None, varvl)
        else:
          varval = var[:,iy,ix]
      
        varst = statsVal(varval)
        var2st = stats2Val(timesv, varval, powerv)
        varstats[iy,ix,0] = varst.minv
        varstats[iy,ix,1] = varst.maxv
        varstats[iy,ix,2] = varst.meanv
        varstats[iy,ix,3] = varst.mean2v
        varstats[iy,ix,4] = varst.stdv
        varquant[iy,ix,:] = varst.quantilesv
        varlregress[iy,ix,:] = var2st.linRegress
        varpregress[iy,ix,:] = var2st.polRegress
        varpregressres[iy,ix] = var2st.polRegressRes
        varpregresssingval[iy,ix,:] = var2st.polRegressSingVal
  elif varshape == 4:
    dz=vardims[1]
    dy=vardims[2]
    dx=vardims[3]
    varstats = np.ones((dz, dy, dx, 5), dtype=float)
    varstats = varstats*varfil
    varquant = np.ones((dz, dy, dx, 21), dtype=float)
    varquant = varquant*varfil
    varlregress = np.ones((dz, dy, dx, 5), dtype=float)
    varlregress = varlregress*varfil
    varpregress = np.ones((dz, dy, dx, powerv+1), dtype=float)
    varpregress = varpregress*varfil
    varpregressres = np.ones((dz,dy,dx), dtype=float)
    varpregressres = varpregressres*varfil
    varpregresssingval = varpregress
    varpregresssingval = varpregresssingval*varfil
    vardnames.append(vardimns[1])
    vardnames.append(vardimns[2])
    vardnames.append(vardimns[3])
    for iz in range(dz):
      for iy in range(dy):
        for ix in range(dx):
          if not varinf.FillValue is None:
            varvl = var[:,iz,iy,ix]
            varval = np.where(varvl == varinf.FillValue, None, varvl)
          else:
            varval = var[:,iz,iy,ix]
      
          varst = statsVal(varval)
          var2st = stats2Val(timesv, varval, powerv)
          varstats[iz,iy,ix,0] = varst.minv
          varstats[iz,iy,ix,1] = varst.maxv
          varstats[iz,iy,ix,2] = varst.meanv
          varstats[iz,iy,ix,3] = varst.mean2v
          varstats[iz,iy,ix,4] = varst.stdv
          varquant[iz,iy,ix,:] = varst.quantilesv
          varlregress[iz,iy,ix,:] = var2st.linRegress
          varpregress[iz,iy,ix,:] = var2st.polRegress
          varpregressres[iz,iy,ix] = var2st.polRegressRes
          varpregresssingval[iz,iy,ix,:] = var2st.polRegressSingVal
  elif varshape == 5:
    dn=vardims[1]
    dz=vardims[2]
    dy=vardims[3]
    dx=vardims[4]
    varstats = np.ones((dn, dz, dy, dx, 5), dtype=float)
    varstats = varstats*varfil
    varquant = np.ones((dn, dz, dy, dx, 21), dtype=float)
    varquant = varquant*varfil
    varlregress = np.ones((dn, dz, dy, dx, 5), dtype=float)
    varlregress = varlregress*varfil
    varpregress = np.ones((dn, dz, dy, dx, powerv+1), dtype=float)
    varpregress = varpregress*varfil
    varpregressres = np.ones((dn,dz,dy,dx), dtype=float)
    varpregressres = varpregressres*varfil
    varpregresssingval = varpregress
    varpregresssingval = varpregresssingval*varfil
    vardnames.append(vardimns[1])
    vardnames.append(vardimns[2])
    vardnames.append(vardimns[3])
    vardnames.append(vardimns[4])
    for iN in range(dn):
      for iz in range(dz):
        for iy in range(dy):
          for ix in range(dx):
            if not varinf.FillValue is None:
              varvl = var[:,iN,iz,iy,ix]
              varval = np.where(varvl == varinf.FillValue, None, varvl)
            else:
              varval = var[:,iN,iz,iy,ix]
      
            varst = statsVal(varval)
            var2st = stats2Val(timesv, varval, powerv)
            varstats[iN,iz,iy,ix,0] = varst.minv
            varstats[iN,iz,iy,ix,1] = varst.maxv
            varstats[iN,iz,iy,ix,2] = varst.meanv
            varstats[iN,iz,iy,ix,3] = varst.mean2v
            varstats[iN,iz,iy,ix,4] = varst.stdv
            varquant[iN,iz,iy,ix,:] = varst.quantilesv
            varlregress[iN,iz,iy,ix,:] = var2st.linRegress
            varpregress[iN,iz,iy,ix,:] = var2st.polRegress
            varpregressres[iN,iz,iy,ix] = var2st.polRegressRes
            varpregresssingval[iN,iz,iy,ix,:] = var2st.polRegressSingVal
  else:
     print errormsg
     print '    timemean: ', varshape, ' shape of matrix not prepared !'
     ncf.close()
     quit(-1)

  vardimnames = tuple(vardnames)
#  print '  shape of desired values: ', vardesiredvalues.shape
# Creation of file
##

  ncfo = NetCDFFile( ofile, 'w')

  vartype = var.dtype
  varattr = var.ncattrs()

# Checking dimensions
##
  newdims = ncfo.dimensions
  Nvardnames = len(vardimnames)
  for idim in range(Nvardnames):
      rdim = vardimnames[idim] 
      if not searchInlist(newdims, rdim):
          if not rdim == 'time':
##              print '      timemean: Adding dimension ' + rdim
              ncfo.sync()
              ncf.close()
              ncfo.close()
              fdimadd(ncfile + ':' + rdim, ofile)
              ncf = NetCDFFile(ncfile,'r')
              ncfo = NetCDFFile(ofile,'a')
          else:
              print '      timemean: No time dimension!'

# Checking fill value
## 
  if searchInlist(varattr, '_FillValue'):
      varfil = var._FillValue
  else:
      varfil = False

  Nstats = len(statsn)
  for ist in range(Nstats):
    newvardnames = []
    if statsn[ist] == 'quant':
      print ist, statsn[ist]##, ': ', varquant[int(dimt/2),10]
      newdim = ncfo.createDimension('quant', 21)
      newvardnames = list(vardnames)
      newvardnames.append('quant')
      newvar = ncfo.createVariable(varn + statsn[ist], vartype, tuple(newvardnames), fill_value=varfil)
      newvar[:] = varquant
    elif statsn[ist] == 'linregress':
      print ist, statsn[ist]##, ': ', varquant[int(dimt/2),10]
      newdim = ncfo.createDimension('lregress', 5)
      newvar = ncfo.createVariable('lregressn', str, ('lregress',))
      newvar[0] = 'slope'
      newvar[1] = 'intercept'
      newvar[2] = 'r_value'
      newvar[3] = 'p_value'
      newvar[4] = 'std_err'
      newvardnames = list(vardnames)
      newvardnames.append('lregress')
      newvar = ncfo.createVariable(varn + statsn[ist], 'f4', tuple(newvardnames), fill_value=varfil)
      newvar[:] = varlregress*1.

    elif statsn[ist] == 'polregress':
      print ist, statsn[ist]##, ': ', varquant[int(dimt/2),10]
      newdim = ncfo.createDimension('pregress', powerv+1)
      newvar = ncfo.createVariable('pregressn', str, ('pregress',))
      for ip in range(powerv+1):
        newvar[ip]='coefficient**' + str(powerv-ip)
      newvardnames = list(vardnames)
      newvardnames.append('pregress')
      newvar = ncfo.createVariable(varn + statsn[ist], 'f4', tuple(newvardnames), fill_value=varfil)
      newvar[:] = varpregress*1.
      newvar.setncattr('power',powerv)
      newvar.setncattr('values','Polynomial coefficients, highest power first')
      newvar = ncfo.createVariable(varn + statsn[ist] + '_Residual', vartype, tuple(vardnames), fill_value=varfil)
      newvar[:] = varpregressres
      newvar.setncattr('power',powerv)
      newvar.setncattr('values','Polynomial residuals')
      newvar = ncfo.createVariable(varn + statsn[ist] + '_VandermondeSingularVector', vartype, tuple(newvardnames), fill_value=varfil)
      newvar[:] = varpregresssingval
      newvar.setncattr('power',powerv)
      newvar.setncattr('values','Polynomial coefficients, highest power first')
    else:
      print ist, statsn[ist]##, ': ', varstats[int(dimt/2),ist]
      if statsn[ist] == 'mean' or statsn[ist] == 'stdv' or statsn[ist] == 'mean2' or statsn[ist] == 'polregress_Residual' \
        or statsn[ist] == 'polregress_VandermondeSingularVector' and searchInlist(varattr, 'scale_factor'):
        newvar = ncfo.createVariable(varn + statsn[ist], 'f4', vardimnames, fill_value=varfil)
        if varshape == 2:
          newvar[:] = varstats[:,ist]*1.
        elif varshape == 3:
          newvar[:] = varstats[:,:,ist]*1.
        elif varshape == 4:
          newvar[:] = varstats[:,:,:,ist]*1.
        elif varshape == 5:
          newvar[:] = varstats[:,:,:,:,ist]*1.
      else:
        newvar = ncfo.createVariable(varn + statsn[ist], vartype, vardimnames, fill_value=varfil)
        if varshape == 2:
          newvar[:] = varstats[:,ist]
        elif varshape == 3:
          newvar[:] = varstats[:,:,ist]
        elif varshape == 4:
          newvar[:] = varstats[:,:,:,ist]
        elif varshape == 5:
          newvar[:] = varstats[:,:,:,:,ist]

    newvar = ncfo.variables[varn + statsn[ist]]
    for attr in varattr:
         newvarattrs = newvar.ncattrs()
         attrv = var.getncattr(attr)
         if not searchInlist(newvarattrs, attr):
              if attr == 'scale_factor' or attr == 'add_offset' or attr == 'valid_range' \
                or attr == 'unpacked_valid_range' or attr == 'actual_range' :
                  if not statsn[ist] == 'mean' and not statsn[ist] == 'stdv' and not statsn[ist] == 'mean2' and not \
                    statsn[ist] == 'linregress' and not statsn[ist] == 'polregress' \
                    and not statsn[ist] == 'polregress_Residual' and not \
                    statsn[ist] == 'polregress_VandermondeSingularVector':
                      newvar.setncattr(attr, attrv)
              else:
                  newvar.setncattr(attr, attrv)
    
    newvar.setncattr('cell_methods', 'time ' + statslongn[ist] + ' all period in file ' + str(dimt) + ' time steps')

  ncfo.sync()
  ncfo.close()

  fvaradd(ncfile + ':lon', ofile)
  fvaradd(ncfile + ':lat', ofile)

  fgaddattr(ncfile, ofile)

  print '    timemean: File "' + ofile + '" as time mean of "' + varn + '" has been created'
