#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys,os
from openmmplumed import PlumedForce 



def load_plumed_script(conf:str):
    """
        load plumed script for openmm simulation
        default units:
        - length: nm
        - energy: kJ/mol
        -time: ps
    """
    script = open(conf['plumed_conf']).read()
    print(script)
    plumed_force = PlumedForce(script=script)
    print('successfully loaded plumed script')
    #plumed_force.setTemperature(conf['temperature']*unit.kelvin)
    return plumed_force

