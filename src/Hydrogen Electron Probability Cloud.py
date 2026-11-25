"""Experimental orbital visualizer; theory context is in ../quantum_mechanics_notes_hydrogen_atom_3d_box.md.

See ../PROJECT_REVIEW.md for the current model limitations and output details.
"""

import numpy as np
import math
import time
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits import mplot3d
from matplotlib import cm
from matplotlib import animation
from itertools import product

#Constants

a_0 = 0.529e-11 #Bohr Radius

#Functions 

def probability(x,y,z,a_0,nlm):
    
    # * math.sin(p) to display 1/4th of the sphere along XZ plane
    # * math.cos(t) to display 1/4th of the sphere along YZ plane
    
    if nlm == '100':
        r = math.sqrt((x**2)+(y**2)+(z**2))
        
        numen = np.abs(math.exp((-2*r)/a_0))
        den = (math.pi*(a_0**3))
        P = (numen/den) * (4*math.pi*math.pow(r,2))
        P = -P
        
    if nlm == '200':
        r = math.sqrt((x**2)+(y**2)+(z**2))
        
        numen = np.abs((math.exp((-r)/a_0)) * math.pow(((2*a_0) - r),2))
        den = np.abs(32*math.pi*math.pow(a_0,5))
        P = (numen/den) * (4*math.pi*math.pow(r,2))
        P = np.sqrt(P)
        P = -P
        
    if nlm == '210':
        eps = 1.0e-10
        r = math.sqrt((x**2)+(y**2)+(z**2))
        t = math.acos(z/(eps + r))
        p = math.atan(y/(x + eps))
        
        numen = np.abs((math.exp((-r)/a_0)) * math.pow(r,2) * math.pow(math.cos(t),2))
        den = np.abs(16*math.pi*math.pow(a_0,5))
        P = (numen/den) * (4*math.pi*math.pow(r,2))
        P = np.sqrt(P)
        
    if nlm == '211':
        eps = 1.0e-10
        r = math.sqrt((x**2)+(y**2)+(z**2))
        t = math.acos(z/(eps + r))
        p = math.atan(y/(x + eps))
        
        numen = np.abs((math.exp((-r)/a_0)) * math.pow(r,2) * math.pow(math.sin(p),2))
        den = np.abs(64*math.pi*math.pow(a_0,5))
        P = (numen/den) * (4*math.pi*math.pow(r,2))
        P = np.sqrt(P)
        
    if nlm == '300':
        r = math.sqrt((x**2)+(y**2)+(z**2))
        
        numen = np.abs((math.exp((-2*r)/3*a_0)) * math.pow(((27*math.pow(a_0,2)) - (18*a_0*r) + (2*math.pow(r,2))),2))
        den = np.abs(6561*math.pi*math.pow(a_0,7))
        P = (numen/den) * (4*math.pi*math.pow(r,2))
        P = np.sqrt(np.sqrt(np.sqrt(np.sqrt(P))))
        P = -P
        
    if nlm == '310':
        eps = 1.0e-9
        r = math.sqrt((x**2)+(y**2)+(z**2))
        t = math.acos(z/(eps + r))
        p = math.atan(y/(x + eps))
        
        numen = np.abs((2*math.pow(r,2)) * (math.exp((-2*r)/3*a_0)) * math.pow((6-(r/a_0)),2) * math.pow(math.cos(t),2))
        den = np.abs(6561*math.pi*math.pow(a_0,5))
        P = (numen/den) * (4*math.pi*math.pow(r,2))
        P = np.sqrt(np.sqrt(np.sqrt(P)))
        P = P
    
    if nlm == '311':
        eps = 1.0e-9
        r = math.sqrt((x**2)+(y**2)+(z**2))
        t = math.acos(z/(eps + r))
        p = math.atan(y/(x + eps))
        
        numen = np.abs((2*math.pow(r,2)) * (math.exp((-2*r)/3*a_0)) * math.pow((6-(r/a_0)),2) * math.pow(math.sin(p),2))
        den = np.abs(6561*math.pi*math.pow(a_0,5))
        P = (numen/den) * (4*math.pi*math.pow(r,2))
        P = np.sqrt(np.sqrt(np.sqrt(P)))
        P = P
        
    if nlm == '320':
        eps = 1.0e-9
        r = math.sqrt((x**2)+(y**2)+(z**2))
        t = math.acos(z/(eps + r))
        p = math.atan(y/(x + eps))
        
        numen = np.abs((math.pow(r,4)) * (math.exp((-2*r)/(3*a_0))) * math.pow(((3*math.pow(math.cos(t),2))-1),2))
        den = np.abs(39366*math.pi*math.pow(a_0,7))
        P = (numen/den) * (4*math.pi*math.pow(r,2))
        P = np.sqrt(np.sqrt(np.sqrt(P)))
        P = P
        
    if nlm == '321':
        eps = 1.0e-9
        r = math.sqrt((x**2)+(y**2)+(z**2))
        t = math.acos(z/(eps + r))
        p = math.atan(y/(x + eps))
        
        numen = np.abs((math.pow(r,4)) * (math.exp((-2*r)/(3*a_0))) * (math.sin(p)**2) * (math.cos(t)**2) * np.real(np.exp(1.j * 2 * t)))
        den = np.abs(6561*math.pi*math.pow(a_0,7))
        P = (numen/den) * (4*math.pi*math.pow(r,2))
        P = np.sqrt(np.sqrt(np.sqrt(P)))
        P = P
    
    if nlm == '322':
        eps = 1.0e-9
        r = math.sqrt((x**2)+(y**2)+(z**2))
        t = math.acos(z/(eps + r))
        p = math.atan(y/(x + eps))
        
        numen = np.abs((math.pow(r,4)) * (math.exp((-2*r)/(3*a_0))) * (math.sin(p)**4) * np.real(np.exp(1.j * 4 * t)))
        den = np.abs(26244*math.pi*math.pow(a_0,7))
        P = (numen/den) * (4*math.pi*math.pow(r,2))
        P = np.sqrt(np.sqrt(np.sqrt(P)))
        P = P
        
    return P

def get_probability_values(x,y,z,a_0,nlm):
    length = len(x)
    L = []
    
    for i in range(length):
        L.append(probability(x[i],y[i],z[i],a_0,nlm))
    
    return np.array(L)

def apply_radius_constraint(coord,n):
    L = []
    R = a_0 * math.pow(n,2)
    for i in coord:
        r = np.sqrt(i[0]**2 + i[1]**2 + i[2]**2)
        if r < R:
            L.append(i)
    return np.array(L)

def set_axes_equal(ax):
    
    #Make axes of 3D plot have equal scale so that spheres appear as spheres,
    #cubes as cubes, etc.

    #Input
    #  ax: a matplotlib axis, e.g., as output from plt.gca().

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def design_orbit(orbit_shape,AzimAngle):
    
    #Bohr radius
    a_0 = 0.529e-11
    
    #Settings
    if orbit_shape == '100':
        f = 2
        resolution = 180
    elif orbit_shape == '200':
        f = 10
        resolution = 200
    elif orbit_shape == '210':
        f = 10
        resolution = 200
    elif orbit_shape == '211':
        f = 10
        resolution = 200
    elif orbit_shape == '300':
        f = 15
        resolution = 160
    elif orbit_shape == '310':
        f = 15
        resolution = 160
    elif orbit_shape == '311':
        f = 15
        resolution = 160
    elif orbit_shape == '320':
        f = 15
        resolution = 160
    elif orbit_shape == '321':
        f = 15
        resolution = 160
    elif orbit_shape == '322':
        f = 15
        resolution = 160
    else:
        print("Invalid Shape Entered")
        
    nx, ny , nz = resolution, resolution, resolution
    scale = f*a_0
    
    #Creating coordinates
    x = np.linspace(-scale,scale, nx)
    y = np.linspace(-scale,scale, ny)
    z = np.linspace(-scale,scale, nz)

    #Creating a 3D cartesian product (3D Mesh) array of dimensions -> (n_x * n_y * n_z) x 3
    res = list(product(x, y, z))
    res = np.array(res)
    n = int(orbit_shape[0])
    res = apply_radius_constraint(res,n)

    #Assigning respective columns to each axis
    x = res[0:(res.shape[0]-1)//2,0]
    y = res[0:(res.shape[0]-1)//2,1]
    z = res[0:(res.shape[0]-1)//2,2]

    #x = res[:,0]
    #y = res[:,1]
    #z = res[:,2]

    
    #Fetching respective probaility value for each coordinate in the 3D Mesh
    probability_values = get_probability_values(x,y,z,a_0,orbit_shape)

    #Plotting the figure
    fig = plt.figure(facecolor="Black", figsize=(6,6))
    ax = fig.add_subplot(projection="3d")

    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    ax.set_box_aspect([1.0, 1.0, 1.0])
    ax.view_init(elev=0.5, azim=-AzimAngle)
    
    # change alpha values of 'rocket' colormap to make low probability regions transparent
    ncolors = 256
    color_array = plt.get_cmap('rocket')(range(ncolors))
    color_array[:,-1] = np.linspace(0.0,0.08,ncolors)
    map_object = LinearSegmentedColormap.from_list(name='rocket_alpha',colors=color_array)
    plt.register_cmap(cmap=map_object)

    im = ax.scatter3D(x, y, z, c=probability_values, marker='.', cmap='rocket_alpha') #sns.color_palette('rocket', as_cmap=True)) #plasma
    
    set_axes_equal(ax)
    ax.set_facecolor("Black")
    ax.set_xlabel('$x-axis$', fontsize=10, rotation=0)
    ax.set_ylabel('$y-axis$', fontsize=10, rotation=0)
    ax.set_zlabel('$z-axis$', fontsize=10, rotation=0)
    
    # Get rid of colored axes planes
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('#000000')
    ax.yaxis.pane.set_edgecolor('#000000')
    ax.zaxis.pane.set_edgecolor('#000000')
    ax.grid(False)
    
    FD = {'fontsize': 10,'color' : 'White'}
    SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    plt.title("Electron Probability Density for Ψ{}".format(orbit_shape.translate(SUB)), fontdict=FD)
    #plt.colorbar(im, label='$|Ψ|^2$', orientation='vertical', ticks=[0,1])
    plt.savefig("Psi{}Images\HydrogenImage{}.png".format(orbit_shape,np.abs(AzimAngle)))
    #plt.show()
    
if __name__ == '__main__':
    
    for angle in range(0,360,1):
        design_orbit('321',angle)
        #design_orbit('321',angle)
    
    #design_orbit('100')
    #design_orbit('200')
    #design_orbit('210')
    #design_orbit('211')
    #design_orbit('300',0.5)
    #design_orbit('310')
    #design_orbit('311')
    #design_orbit('320')
    #design_orbit('321')
    #design_orbit('322')
    
    
    
    