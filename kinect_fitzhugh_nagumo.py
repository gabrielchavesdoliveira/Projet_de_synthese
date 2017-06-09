#Simulation of the Fitzhguh-Nagumo reaction diffusion system





import freenect
import cv2
import sys
import random
import numpy as np
import time
import scipy
from freenect import sync_get_depth as get_depth
import fitzhugh_nagumo



print(" Press : ")
print("   s : start/pause")
print("   i : reinitialize the concentrations")
print("   q : quit")
print("   c : erase the reactant v in a randomly chosen box patch")
print("   m : mask the reactant with a randomly generated mask")
print("   p : save the current u potential")
print("   f : toggle fullscreen/normal screen")


try:
    fullscreen_flag = cv2.WINDOW_FULLSCREEN
    normal_flag = cv2.WINDOW_NORMAL
except:
    fullscreen_flag = cv2.cv.CV_WINDOW_FULLSCREEN
    normal_flag = cv2.cv.CV_WINDOW_NORMAL

cv2.namedWindow('u', cv2.WND_PROP_FULLSCREEN)
cv2.setWindowProperty("u", cv2.WND_PROP_FULLSCREEN, cv2.WINDOW_NORMAL)


key = 0
run = False

#
d = 3.0
width = 200
height = 100
dt = 1
pattern = 'propagation'
display_scaling_factor = 4
# The frustum for the kinect depth
zmin = 1
zmax = 4




model = fitzhugh_nagumo.SpectralModel(pattern, width=width, height=height, d=d, dt=dt, mode='ETDFD')


model.init()





def make_effect(u_orig, scale):
    res_height, res_width = scale * u_orig.shape[0], scale * u_orig.shape[1]
    s_kernel = 11
    kernel = np.ones((s_kernel,s_kernel), dtype=np.float)
    # Light coming from top left
    kernel[:int(2./3 * s_kernel),:int(2./3 * s_kernel)] = -1
    # Light coming from left
    #kernel[:int(2./3 * s_kernel),:] = -1
    effect = scipy.signal.convolve2d(2. * (u_orig - 0.5), kernel, mode='same')
    effect /= 30. # HAND TUNED SCALING of the effect ... might need to be adapted if changing s_kernel
    effect = 0.8*effect
    effect[effect >= 1.0] = 1.0
    effect[effect <= 0.0] = 0.0
    effect_hires = cv2.resize(effect, (res_width, res_height), interpolation=cv2.INTER_CUBIC)

    u_hires = cv2.resize(u_orig, (res_width, res_height),interpolation=cv2.INTER_CUBIC)
    # Blur the image to get the shading
    u_blur = scipy.ndimage.filters.uniform_filter(u_hires, size=5)
    # Shift the shadding down right
    u_blur = np.lib.pad(u_blur, ((2,0),(2,0)), 'constant', constant_values=1)[:-2,:-2]
    
    dst = 0.6 * u_hires + 0.4 * effect_hires
    dst[u_hires >= 0.99] = u_blur[u_hires >= 0.99]
    

    dst[dst<2] = 1 - dst[dst<2]
            
    return dst

def insert_text(img, text):
    global height
    img[(img.shape[0]-40):,:] = 1
    cv2.putText(img, text, (20, img.shape[0]-10), cv2.FONT_HERSHEY_SIMPLEX, 1, 0, thickness=2)

def insert_depth(depth_img, img):
    tgt_size = (2*64, 2*48)
    small_depth = cv2.resize(depth_img, tgt_size)
    #print(small_depth.shape)
    img[(img.shape[0]-small_depth.shape[0]):, (img.shape[1]-small_depth.shape[1]):] = small_depth[:,:, 0]
    return



u = np.zeros((height, width))
depth = np.ones((height, width))
depth_img = np.zeros((2,2,3), dtype=np.float)
can_mask = False
epoch = 0
N = (height*width)**0.5 
t0 = time.time()
frame_id = 0
#format_video=cv2.cv.CV_FOURCC(*'XVID')

#nom_video="video1"

#video=cv2.VideoWriter(nom_video, format_video, 23, (200,100))

while key != ord('q'):
    # As soon as we get a minimum reactant, we start
    # to take into account the kinect
    # Otherwise, the whole activity vanishes
    if((model.get_ut().mean() <= 0.9) and not can_mask):
        can_mask = True
        print("Masking begins")
    if(run):
        epoch += 1
        #print(epoch)
        if(epoch == 100):
            t1 = time.time()
	    print("FPS: %f fps" % (100 / (t1 - t0)))
            print("Time = " + str(t1))
	    t0 = t1
            epoch = 0
        if(can_mask and (epoch % 2 == 0)):
            (depth,_) = get_depth(format=freenect.DEPTH_MM)
            # Restrict to box in [zmin; zmax]
            # depth is scaled in meters, and the horizontal axis is flipped
            # what is in [zmin, zmax] is rescaled to [1, 0] , the rest set to 0
            depth = (zmax - depth[:,::-1]*1e-3)/(zmax - zmin)
            depth[depth < 0] = 0
            depth[depth > 1] = 0

            depth_img = (np.dstack((depth, depth, depth)).astype(np.float))
            #cv2.resize(depth_img, (width, height))
            #cv2.imshow('Depth', depth_img)
            
            depth = cv2.resize(depth.astype(np.float), (width, height))
	    #depth = 1. - cv2.resize(depth, (N, N))
            #print(depth.min(), depth.max(), depth.mean())
            #depth = depth * 0.85 / depth.mean()
            #mask = 0.75 + 0.25 * depth
        model.mask_reactant(depth)
        model.step()
        u[:,:] = model.get_ut()

    u_img = make_effect(u, display_scaling_factor)
    insert_text(u_img, "Fitzhugh-Nagumo Reaction Diffusion")
    insert_depth(depth_img, u_img)

    cv2.imshow('u', u_img)
    
    #u_img2 = np.dstack(255*(u_img,u_img,u_img)).astype(np.float)
    #u_img2[u_img>255] = 255
    #u_img2[u_img<0] = 0
    

    key = cv2.waitKey(1) & 0xFF
    #video.write(u_img2) 
    
    if(key == ord('c')):
        c = (random.randint(0, N-1), random.randint(0, N-1))
        model.erase_reactant(c , N/8)
    elif(key == ord('m')):
        mask = 0.75 + 0.25*np.random.random((N, N))
        model.mask_reactant(mask)
    elif key == ord('s'):
        run = not run
        print("Running ? : " + str(run))
    elif key == ord('i'):
        model.init()
    elif key == ord('p'):
        print("Saving u-%05d.png" % frame_id)
        cv2.imwrite("u-%05d.png" % frame_id, (np.minimum(255*u_img, 255)).astype(np.uint8))
        frame_id += 1
    elif key == ord('f'):
        screenmode = cv2.getWindowProperty("u", cv2.WND_PROP_FULLSCREEN)
        if(screenmode == normal_flag):
            cv2.setWindowProperty("u", cv2.WND_PROP_FULLSCREEN, fullscreen_flag)
        else:
            cv2.setWindowProperty("u", cv2.WND_PROP_FULLSCREEN, normal_flag)


#video.release()
#cv2.destroyAllWindows()
