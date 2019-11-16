# -*- coding: utf-8 -*
"""
Created on Sun Oct 20 12:26:08 2019

@author: fghjk
"""

from uwimg import *
#im = load_image("data/dog.jpg")
#f = make_box_filter(7)
#blur = convolve_image(im, f, 1)
#thumb = nn_resize(blur, blur.w//7, blur.h//7)
#save_image(thumb, "dogthumb")
#make_highpass_filter
im = load_image("data/dog.jpg")
f = make_highpass_filter()
blurs = convolve_image(im, f, 0)
clamp_image(blurs)
test=load_image("figs/dog-highpass.png")
blurr=convolve_image(test,f,0)
save_image(blurs, "dog_highpass_filter")
same_image(blurs,test)
#spharp
#im = load_image("data/dog.jpg")
#f = make_sharpen_filter()
#blurs = convolve_image(im, f, 1)
#save_image(blurs, "dog_sharpen")
##eboss
#im = load_image("data/dog.jpg")
#f = make_emboss_filter()
#blurs = convolve_image(im, f, 1)
#test=load_image("figs/dog-emboss.png")
##save_image(blurs, "dog_emboss")
#same_image(blurs,test)
#gauss
#im = load_image("data/dog.jpg")
#f = make_gaussian_filter(2)
#blur = convolve_image(im, f, 1)
#save_image(blur, "dog-gauss2")
##low high
#im = load_image("data/dog.jpg")
#f = make_gaussian_filter(2)
#lfreq = convolve_image(im, f, 1)
#hfreq = im - lfreq
#reconstruct = lfreq + hfreq
#save_image(lfreq, "low-frequency")
#save_image(hfreq, "high-frequency")
#save_image(reconstruct, "reconstruct")
##sobel
#im = load_image("data/dog.jpg")
#res = sobel_image(im)
#mag = res[0]
#feature_normalize(mag)
#save_image(mag, "magnitude")
#colorsobel
im = load_image("data/dog.jpg")
res = colorize_sobel(im)


save_image(res, "magnitudec")