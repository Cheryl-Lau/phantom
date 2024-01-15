# -*- coding: utf-8 -*-
"""
For deleting dumpfiles 
"""

import os 

pwd = os.getcwd()
for file in os.listdir(pwd):
    if file.startswith("sphere_"):
        print(file)
        
        ends = ( "1","2","3","4","5","6","7","8","9")

        if file.endswith(ends):
            os.remove(file)