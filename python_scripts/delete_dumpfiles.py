# -*- coding: utf-8 -*-
"""
For deleting dumpfiles 
"""

import os 
import re

dump_prefix = "sphere_"
idump_lastdumptocheck = 20

print('dump prefix: ',dump_prefix)
print('last idump to delete: ',idump_lastdumptocheck)

proceed = input('Proceed [y/n]? ')
if (proceed.strip() != 'y'):
    quit() 

pwd = os.getcwd()
for file in os.listdir(pwd):
    if file.startswith(dump_prefix.strip()):
        idump = re.findall(r'\d+',file.replace(dump_prefix, ''))
        idump = int(idump[0])

        if (idump < idump_lastdumptocheck):

            ends = ( "1","2","3","4","5","6","7","8","9")

            if file.endswith(ends):
                print('Deleting ',file)
                os.remove(file)
