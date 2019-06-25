# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 16:56:29 2019

@author: sarung
"""

from skimage.segmentation import slic
from skimage.segmentation import mark_boundaries
import skimage
import matplotlib.pyplot as plt

b = skimage.data.camera()
segments = slic(b, n_segments=1000, compactness=0.1, sigma=1)
plt.imshow(mark_boundaries(b, segments))
plt.show()

#%%
from skimage.io import imread, imsave
import numpy as np
#fdir = 'C:\\Users\\sarung\\Documents\\LipingsData\\ZstackL1_3_2grayscale'
#t = 220
#imgs = []
#for z in range(9):
#    im = imread(f'{fdir}\\larva3_2_z{z:01d}_t{t:03d}.tif')
#    im -= im.min()
#    im = im / im.max()
#    imgs.append(im)
#V = np.stack(imgs)
#plt.imshow(V.max(axis=0), cmap='gray')
#imsave(f'larva3_2_t{t:03d}.png', V.max(axis=0))


#%%
from skimage.transform import resize
from skimage import color

s = (100,100)
im = imread('larva3_2_t218_s1.png')
im = resize(im, s)
#plt.imshow(im)
segments = slic(im, n_segments=50, compactness=100)
#imsave(f'input.png', im)

im = color.rgb2gray(im)
imb = mark_boundaries(im, segments, color=(237/255, 190/255, 0))
#imsave(f'frame_img.png', imb)
#%%
imslic = np.zeros(s)
for i in range(segments.max()):
    mask = segments == i
    imslic[mask] = np.mean(im[mask])

imb = mark_boundaries(imslic, segments, color=(237/255, 190/255, 0))
#imsave(f'frame_img2.png', imb)

imslic = imslic > .5
#imslic -= imslic.min()
#imslic = imslic / imslic.max()

imb = mark_boundaries(imslic, segments, color=(237/255, 190/255, 0))
plt.imshow(imslic)
plt.show()
imsave(f'template_slic.png', imb)
imsave(f'template_img.png', imslic*255)

#


#%%
#from skimage.transform import resize
#from skimage import color
#
#s = (100,100)
#im = imread('larva3_2_t231_s1.png')
#im = resize(im, s)
##plt.imshow(im)
#segments = slic(im, n_segments=50, compactness=100)
#
#im = color.rgb2gray(im)
#imb = mark_boundaries(im, segments, color=(237/255, 190/255, 0))
#
#plt.imshow(imb)
#imsave(f'frame_img_t231.png', imb)