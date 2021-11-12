import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import cv2 as cv
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
import matplotlib.image as mpimg

for alpha in range(0,11):
    os.system("./sot -n 2000 -p 2000 --step 32 -d 3 -f1 ./pig.obj -f2 ./dog.obj -a "+str(alpha)+" -o ./pointsets/"+str(alpha)+".dat")
    os.system("python3 ../../../utils/render_pointset.py ./pointsets/"+str(alpha)+".dat ./results_morphing/"+str(alpha)+".png")