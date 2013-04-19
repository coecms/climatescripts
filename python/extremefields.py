#!/usr/bin/env python
"""
file:   extremefields.py
author: Scott Wales (scott.wales@unimelb.edu.au)

Checks a data file for extreme values in the component fields by measuring
the number of standard deviations the min and max of each field is from the
mean.

Uses cdat-lite to open the file, it will work for netcdf files but may have
issues with UM formatted files. If you do have trouble convert the UM file to
netcdf using xconv.

To run at NCI:

    module use /g/data/access/modules
    module load cdat-lite
    extremefields.py FILE

================================================================================

Copyright 2013 ARC Centre of Excellence for Climate System Science

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import numpy
import cdms2
import math
import argparse

parser = argparse.ArgumentParser(description="Check extreme values in a climate data file by printing the number of standard deviations each field's min and max is from the field mean")
parser.add_argument('filename', metavar='FILE', nargs=1, help="File to analyse")
args = parser.parse_args()
print args.filename

file = cdms2.open(args.filename[0])

extremes = {}

# Create a tuple for each variable in the file containing the number of
# standard deviations the min and max are from the mean
for var in file.variables:
    field = file[var]
    mean = numpy.mean(field)
    std = numpy.std(field)

    extremes[var] = ((numpy.min(field)-mean)/std, (numpy.max(field)-mean)/std)

# Sort the values by the sum of squares, ignoring nan's
sortedExtremes = sorted(extremes.items(), 
        key=lambda x: float('-inf') if math.isnan(x[1][0]) else x[1][0]**2+x[1][1]**2)

# Print the values
for value in sortedExtremes:
    print "% 60s\t% 8.2e\t% 8.2e"%(file[value[0]].attributes['long_name'],value[1][0],value[1][1])
