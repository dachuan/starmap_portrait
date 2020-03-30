#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from math import modf

import astropy.time
from astropy.coordinates import EarthLocation,SkyCoord,get_sun,get_moon,get_body
import astropy.units as u

from datetime import datetime, timedelta
import time

from geopy.geocoders import Nominatim


#
# handle observation time
#

#obs_time = astropy.time.Time( datetime.utcnow() )

custom_time = datetime(2020, 3, 30, 0, 0, 0, 00000)
utc_offset = (time.timezone if (time.localtime().tm_isdst==0) else time.altzone)/60/60*-1
obs_time = astropy.time.Time(custom_time - timedelta(hours=modf(utc_offset)[1],
                                                     minutes=modf(utc_offset)[0]))
#
# handle observation location
#

# set manually
#latitude = 31 #31
#longitude = 121 #121
#altitude = 1
#obs_loc = EarthLocation( lat=latitude*u.deg,
#                             lon=longitude*u.deg,
#                             height=altitude*u.m )

geolocator = Nominatim(user_agent='skymap')
location = geolocator.geocode("Leiden", timeout=5)
assert location is not None, 'location is not set'

obs_loc = EarthLocation( lat=location.latitude*u.deg,
                         lon=location.longitude*u.deg,
                         height=location.altitude*u.m )
#
## set main plt para
#

plt.style.use('fivethirtyeight')
figsize = (14,14) #4.8,4.8
dpi = 160 #300

fonts ={
        'cons':5, # 4.5
        'star':4.5, #4.5
        'planet':4,
        }

cnfont = fm.FontProperties(fname='/System/Library/Fonts/PingFang.ttc')
iscn = False # chinese or english name of entities

colors ={
        'sky':'#000000',
        'fov':'#8899CC',
        'cons_line':'#FFFFFF',
        'zodiac':'#FFC704',
        'cons_name':'#FFFFFF',
        'star_name':'#FFFFFF',
        'solar_name':'#FFC704',
        }

alphas ={
        'cont':0.8,
        'star':0.8,
        'solar':0.6,
        'cont_line':0.9,
        'fov':0.4
        }

linewidth = 0.6   # 0.6 constellations lines
mainstar_size_scalar = 2.5 # scale the size of star by mag
star_size_scalar = 1.2 #1 scale the size of star by mag
star_halo_scalar = 40 #20 halo size of star by mag

#
## observation horison is (0,0)(0,180)
## change from alt az to ra dec
#
def obs_horison(obs_time,obs_loc):
    SN_ax_alt = [0, 0]
    SN_ax_az = [0, 180]
    SN_ax = SkyCoord( SN_ax_az, SN_ax_alt, unit='deg', frame='altaz',\
                      obstime=obs_time, location=obs_loc )
    SN_ax = SN_ax.transform_to('icrs')
    # SN_ax[0], SN_ax[1]
    # todo: return varys a lot 
    return SN_ax

#-------------------- important for all transform--------------------# 

global_SN_ax = obs_horison(obs_time,obs_loc) # global

#--------------------------------------------------------------------# 


#
# make plot field
#

fig = plt.figure(dpi=dpi)
fig.canvas.set_window_title('starmap-northern')

# import image as background
img = plt.imread('./res/rhine.jpg') #需要安装pillow
ax_image = fig.add_axes([0,0,1,1],label='img')
ax_image.imshow(img,alpha=1)
ax_image.axis('off')
print("portrait is ready")

# set skymap ax
ax = fig.add_axes([0, 0, 1, 1], polar=True,label='sky')
ax.set_theta_direction(-1)  # anti-clockwise, according to  stars' ra rule
ax.set_xticks([]) #no ticks
ax.set_yticks([]) #no ticks
ax.set_ylim(-45, 45-global_SN_ax[1].dec.degree + 1) # 1 --margin

ax.set_facecolor(colors['sky'])

# set polar ax transparent
ax.patch.set_alpha(0)
print("sky map is prepared")

#fig,ax = prepare_obs_skymap()

#
## curve text for display star, constellation names
#

def curve_txt(theta0,r,s,ax,ditch=3.2,**kw):

    #ditch = 3.2 # 3.6 distance between letters
    d_theta = ditch / (r+45) # -45 is the polar point

    # check theta0, the first letter pos
    
    if theta0 < 0:
        theta0 += 2*np.pi # negtive to positive

    if theta0 > np.pi:
        s = s[::-1] # if upper screen , reverse the string of name

    for i,v in enumerate(s):
        angle = theta0 - i*d_theta

        # upper screen and bottom is different in rotation
        if theta0 < np.pi:
            r_angle = angle * 180 /np.pi - 90
        else:
            r_angle = angle * 180 /np.pi - 270

        ## ax.set_theta_direction(-1)
        r_angle = - r_angle

        ## set text
        ax.text(angle, r, v,
                rotation=r_angle, rotation_mode='anchor',
                va='center',ha='center',
                **kw)



#
## constellations data
#

asterisms= pd.read_csv('./res/58constellations.csv') # 58aster with mag

## rotation axe according to the observation
ax_rotation_angle = global_SN_ax[0].ra.radian

## for check name position
## tmp
#ax_rotation_angle = 0

for index,row in asterisms.iterrows():

    # ra in source data is 24h 
    ras = [float(x)*360/24 for x in row['ra'].replace('[', '').replace(']', '').split(',')]
    decs = [float(x) for x in row['dec'].replace('[', '').replace(']', '').split(',')]

    ras = [np.radians(i) for i in ras]
    decs = [-i + 45 for i in decs]

    ## rotated in new ra
    ## angle - 90(North face) - abservation ratation angle
    ## avoid using axe rotation, set_theta_orieation
    ras = [i-ax_rotation_angle-np.pi/2 for i in ras]

    ## plot constellations lines
    #------------------------------------------------------------#

    if row['zodiac'] == True:
        cons_lc = colors['zodiac']
    else:
        cons_lc = colors['cons_line']

    for n in range(int(len(ras)/2)): #1-2
        ax.plot(ras[n*2:(n+1)*2],decs[n*2:(n+1)*2],
                color=colors['cons_line'], lw=linewidth,
                alpha=alphas['cont_line'],zorder=0)



    ## plot stars in constellations
    #------------------------------------------------------------#
    # set star size by mag
    mags = row['mags'].replace('[','').replace(']','')\
            .replace(' ','').split(',')
    mags = [float(i) for i in mags]

    st_v = 5.  # ________/---------- valve

    for i,mag in enumerate(mags):

        if mag < 0.1:
            ss = st_v - 0.1
        else:
            ss = st_v - mag

        if ss < 0:
            ss = 0

        ax.scatter(ras[i], decs[i],
           s=ss*mainstar_size_scalar, color='white', lw=0, edgecolor='none', 
           alpha=0.8, zorder=1) 


## plot the name of constellations
#------------------------------------------------------------
#
# plot constellations name
## read pos and name from txt
# 
df_pos = pd.read_csv('./res/con_name_pos.txt', sep=':',header=None)
df_pos.columns = ['ra','dec','en','cn']

cn_cont = iscn 

if cn_cont:
    for index,row in df_pos.iterrows():
        curve_txt(-ax_rotation_angle + row['ra']/180*np.pi,row['dec'],row['cn'],ax,
                ditch = 5,
                color=colors['cons_name'],
                FontProperties=cnfont,fontsize=fonts['cons'],
                alpha=alphas['cont'],
                weight='semibold', #'roman'
                zorder=5)
else:
    for index,row in df_pos.iterrows():
        curve_txt(-ax_rotation_angle + row['ra']/180*np.pi,row['dec'],row['en'],ax,
                color=colors['cons_name'],fontsize=fonts['cons'],
                alpha=alphas['cont'],
                weight='semibold', #'roman'
                zorder=5)

#
## plot all sky stars 
#------------------------------------------------------------
#
## data for stars
#
stars = pd.read_csv('./res/hygdata_processed_mag65.csv')
stars_nonvar = stars[pd.isnull(stars['var'])]
stars_nonvar = stars_nonvar[stars_nonvar['color'] != '#000000']

#
# plot star names
#
star_pos = pd.read_csv('./res/star_name_pos.txt', sep=':',header=None)
star_pos.columns = ['ra','dec','en','cn']

cn_star = iscn

if cn_star:
    for index,row in star_pos.iterrows():
        curve_txt(-ax_rotation_angle + row['ra']/180*np.pi,row['dec'],row['cn'],ax,
                ditch = 5,
                color=colors['star_name'],
                FontProperties=cnfont,fontsize=fonts['star'],
                alpha=alphas['star'],
                weight='semibold', #'roman'
                zorder=6)
else:
    for index,row in star_pos.iterrows():
        curve_txt(-ax_rotation_angle + row['ra']/180*np.pi,row['dec'],row['en'],ax,
                color=colors['star_name'],fontsize=fonts['star'],
                alpha=alphas['star'],
                weight='semibold', #'roman'
                zorder=6)

#  
# plot stars, seg by mag
#
mag47 = stars_nonvar[stars_nonvar['mag'] < 4.7]
mag3 = stars_nonvar[stars_nonvar['mag'] < 3]

# stars in mag47 is 879, a little bit slow,
for index, row in mag47.iterrows(): 
        # fill the star
        # zorder is the rendering order z-index-order
        # todo: may use spec color=row['color']
   mg_v = 4.7
   if row['mag'] < 3:
       ss = mg_v - 3
   else:
       ss = mg_v - row['mag']

   ax.scatter(np.radians(row['ra']*360/24)-ax_rotation_angle-np.pi/2, 45 - row['dec'],
      s=ss*star_size_scalar, color='white', lw=0, edgecolor='none', 
      alpha=0.8, zorder=2)

# halo like as light3 star
for index, row in mag3.iterrows(): 
   #alpha=min(1, (mg_v-row['mag'])/mg_v),
   mg_v = 3
   ax.scatter(np.radians(row['ra']*360/24)-ax_rotation_angle-np.pi/2, 45 - row['dec'],
      s=(mg_v - row['mag'])*star_halo_scalar, color='white', lw=0, edgecolor='none', 
      alpha=0.15,  zorder=3)

#
# plot current field-of-view
#

def plot_obs_fov(ax, obs_time, obs_loc,ax_rotation_offset):

    # coordinates of field-of-view-circle' points in (alt,az) frame
    fov_az = np.arange(0, 360+0.1, 1)
    fov_alt = np.zeros(len(fov_az))
    fov = SkyCoord(fov_az, fov_alt, unit='deg', frame='altaz', obstime=obs_time, location=obs_loc)

    # converting this coordinates to (RA,dec) format for plotting them onto plot
    fov = fov.transform_to('icrs')

    # fill the area that we cannot observe now
    shared_ax = fov.ra.radian -np.pi/2 - ax_rotation_offset
    fov_circle = -fov.dec.value+45

    # set outside circle 
    out_circle_d = 90 - global_SN_ax[1].dec.degree +1 #1 --- margin
    outer_circle = len(fov_circle) * [out_circle_d]
    ax.fill_between( shared_ax, fov_circle, outer_circle, where=outer_circle>=fov_circle,
                     facecolor=colors['fov'], alpha=alphas['fov'],zorder=10 )

    print("field-of-view is plotted")

#
# plot sun
#
def plot_sun(ax,obs_time):
    
    sun = get_sun(obs_time)
    sun = SkyCoord(sun.ra, sun.dec, frame='gcrs').transform_to('icrs')
    theta = sun.ra.radian - ax_rotation_angle -np.pi/2
    r = 45-sun.dec.value

    # body
    ax.scatter(theta, r, color='white', s = 40, zorder=7)
    # halo
    ax.scatter(theta, r, color='white', s = 40*25, zorder=7,alpha=0.03)
    ax.scatter(theta, r, color='white', s = 40*20, zorder=7,alpha=0.05)
    ax.scatter(theta, r, color='white', s = 40*9,  zorder=7,alpha=0.05)
    # text
    cn = iscn
    if cn:
        fp = cnfont
        txt = '太 阳'
        curve_txt(theta+np.pi/180*2,r+10,txt,ax,
                  weight='semibold',color=colors['solar_name'],zorder=7,
                  FontProperties=fp,fontsize=6,
                  alpha=alphas['solar'])
    else:
        txt = 'Sun'
        curve_txt(theta+np.pi/180*2,r+10,txt,ax,
                  weight='semibold',color=colors['solar_name'],zorder=7,
                  fontsize=6,
                  alpha=alphas['solar'])

    print('sun is plotted.')
    
#
# plot moon ,pos with phase
#
def moon_phase_angel(obs_time):

    sun = get_sun(obs_time)
    moon = get_moon(obs_time)

    elongation = sun.separation(moon)
    return np.arctan2(sun.distance*np.sin(elongation),
            moon.distance-sun.distance*np.cos(elongation))
    
def moon_illumination(obs_time):
    i = moon_phase_angel(obs_time)
    k = (1 + np.cos(i))/2.0
    return k.value


def plot_moon(ax,obs_time,obs_loc):
    
    moon = get_moon(obs_time, location=obs_loc)
    moon = SkyCoord(moon.ra, moon.dec, frame='gcrs').transform_to('icrs')
    theta = moon.ra.radian - ax_rotation_angle -np.pi/2
    r = 45-moon.dec.value

    # body
    #ax.scatter(theta, r, color='white', s = 25, zorder=7)
    # halo
    ax.scatter(theta, r, color='white', s = 25*25, zorder=7,alpha=0.03)
    ax.scatter(theta, r, color='white', s = 25*20, zorder=7,alpha=0.05)
    ax.scatter(theta, r, color='white', s = 25*9,  zorder=7,alpha=0.05)

    # text
    cn = iscn
    if cn:
        txt = '月亮'
        curve_txt(theta+np.pi/180*2,r+5,txt,ax,
                  weight='semibold',color=colors['solar_name'],zorder=7,
                  FontProperties=cnfont,
                  fontsize=4,
                  alpha=alphas['solar'])
    else:
        txt = 'moon'
        curve_txt(theta+np.pi/180*2,r+5,txt,ax,
                  weight='semibold',color=colors['solar_name'],zorder=7,
                  fontsize=4,
                  alpha=alphas['solar'])

    # moon phase
    # ( --->  ( is a turn
    # ?% percent of this turn
    illumination = moon_illumination(obs_time)
    #print("moon circle per: ", illumination)

    # use moon_font.ttf
    # change fraction to a symbol
    symbol = illumination * 26 # 26 letters from A - Z
    if symbol<0.2 or symbol > 25.8:
        symbol = '1'
    else:
        symbol = chr(ord('a') + int(symbol + 0.5) - 1)
    #print(symbol)

    # plot symbol, 
    prop = fm.FontProperties(fname='./res/moon_phases.ttf',size=8)
    ax.text(theta,r,symbol,
            ha='center', va='center',
            color='white',FontProperties=prop)

    print('moon is plotted.')


#
# plot planets
#
def plot_planets(ax,obs_time,obs_loc):

    planets = ['Mercury',
               'Venus',
               'Mars',
               'Jupiter',
               'Saturn',
               #'Uranus',
               #'Neptune'
               ]

    cn = iscn
    if cn:
        names = ['水星',
                    '金星',
                    '火星',
                    '木星',
                    '土星',
                    #'天王星',
                    #'海王星'
                    ]
    else:
        names = ['Mercury',
                    'Venus',
                    'Mars',
                    'Jupiter',
                    'Saturn',
                    #'Uranus',
                    #'Neptune'
                    ]

    planets_coords = [get_body(planet, obs_time, location=obs_loc) for planet in planets]

    for name,coords in zip(names,planets_coords):

        planet = SkyCoord(coords.ra, coords.dec, frame='gcrs').transform_to('icrs')
        theta = planet.ra.radian - ax_rotation_angle -np.pi/2
        r = 45 - planet.dec.value
    
        # body
        ax.scatter(theta, r, color='white', s = 12, zorder=7)
        # halo
        ax.scatter(theta, r, color='white', s = 12*25, zorder=7,alpha=0.03)
        ax.scatter(theta, r, color='white', s = 12*20, zorder=7,alpha=0.05)
        ax.scatter(theta, r, color='white', s = 12*9,  zorder=7,alpha=0.05)
        # text
        if cn:
            curve_txt(theta+np.pi/180*2,r+4, name, ax,
                    weight='semibold',
                    FontProperties=cnfont,fontsize=fonts['planet'],
                    color=colors['solar_name'],zorder=7,
                    alpha=alphas['solar'])
        else:
            curve_txt(theta+np.pi/180*2,r+4, name, ax,
                    weight='semibold',
                    fontsize=fonts['planet'],
                    color=colors['solar_name'],zorder=7,
                    alpha=alphas['solar'])


    print('planets are plotted.')




# plot all entities
plot_obs_fov(ax, obs_time, obs_loc,ax_rotation_angle)
plot_sun(ax,obs_time)
plot_moon(ax,obs_time,obs_loc)
plot_planets(ax,obs_time,obs_loc)

# tmp
#ax.axis('off')


fig.set_size_inches(4.8,4.8)
fig.savefig('./pro/out.png',dpi=300)
#fig.savefig('./pro/out.png',dpi=300,transparent=True) # transparent background

plt.show()
