#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import os, sys

def harmonic_cos_angle_form(use_pbc, force_group=2):
	"""
	cosine/squared-style harmonic potential:
	U_angle = \sum K_angle (\cos \theta - \cos \theta_0)^2
	"""
	
	angle_expression = 'K * (cos(theta) - cos(theta0))^2'
	angle_force      = mm.CustomAngleForce(angle_expression)
	angle_force.addPerAngleParameter("theta0")
	angle_force.addPerAngleParameter("K")

	

