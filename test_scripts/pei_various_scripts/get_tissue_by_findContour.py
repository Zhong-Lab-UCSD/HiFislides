import sys
import numpy as np
import cv2
import pandas as pd
from skimage import morphology
from scipy.ndimage import binary_fill_holes

infile = sys.argv[1]
coord = sys.argv[2]
myt = sys.argv[3]

img = cv2.imread(infile,cv2.IMREAD_COLOR)
 
mycolnames_0 = ["X","Y"]
tar = pd.read_csv(coord,sep="\t",names=mycolnames_0)

(h, w) = img.shape[:2]

# OpenCV reads color images as BGR instead of RGB
# OR: rgb_image = cv2.cvtColor(img,cv2.COLOR_BGR2RGB)

img_grey = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)

#set a thresh
thresh = int(myt)

#get threshold image
ret,thresh_img = cv2.threshold(img_grey, thresh, 255, cv2.THRESH_BINARY)

#find contours
contours, hierarchy = cv2.findContours(thresh_img, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

sorted_contours = sorted(contours, key=cv2.contourArea, reverse=True)

cnti = 1;
for cnt in sorted_contours:
	k = 0
	for i,rowi in tar.iterrows():
		# rowi["X"]: from left to right.
		# rowi["Y"]: from top to down.
		dist = cv2.pointPolygonTest(cnt,(float(rowi["X"]),float(rowi["Y"])),True)
		if dist > 0:
			k = k + 1
	if k == tar.shape[0]:
		area = cv2.contourArea(cnt)
		perimeter = cv2.arcLength(cnt,True)
		print(cnti,k,area,perimeter,"Hit!\n")
		img2 = img
		opfile = sys.argv[4] + "_" + myt + "_" + str(cnti) + ".png"
		for ii in range(0,h):
			for ij in range(0,w):
				hit = 0
				dist = cv2.pointPolygonTest(cnt,(ij,ii),True)
				if dist > 0:
					# row number first, then column number.
					img2[ii,ij] = img[ii,ij]
					hit = 1
				if hit == 0:
					img2[ii,ij] = (255,255,255)
		cv2.imwrite(opfile,img2)
	cnti = cnti + 1
