# -*- coding: utf-8 -*-
"""
For deleting dumpfiles 
"""

import os 
import re

dump_prefix = "sphere_"
idump_lastdumptocheck = 1000

pwd = os.getcwd()
for file in os.listdir(pwd):
    if file.startswith(dump_prefix.strip()):
        idump = re.findall(r'\d+',file)
        idump = int(idump[0])
        print(file,idump)

        if (idump < idump_lastdumptocheck):

            ends = ( "1","2","3","4","5","6","7","8","9")

            if file.endswith(ends):
                os.remove(file)
