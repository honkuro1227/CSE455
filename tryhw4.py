# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 06:43:14 2019

@author: fghjk
"""
from uwimg import *

a = load_image("data/dog_a.jpg")
b = load_image("data/dog_b.jpg")

flow = optical_flow_images(b, a, 15, 8)
draw_flow(a, flow, 8)
save_image(a, "linesr")
optical_flow_webcam(15,4,8)