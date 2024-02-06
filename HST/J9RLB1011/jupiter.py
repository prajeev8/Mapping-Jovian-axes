#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 12:29:18 2024

@author: pallavirajeev
"""

#Figure 1

# from astropy.io import fits
# import matplotlib.pyplot as plt
# from astropy.wcs import WCS
# # from astropy.coordinates import SkyCoord
# from astropy.coordinates import Angle

# def plot_file(file_name):
#     with fits.open(file_name) as hdu_list:
#         header = hdu_list[0].header + hdu_list[1].header
#         w = WCS(header)
#         fig, ax = plt.subplots()
#         plt.imshow(hdu_list[1].data, origin='lower')

#         # Calculate pixel coordinates based on WCS information
#         ax_pix, ay_pix = w.wcs.crval
#         x_pix, y_pix = w.wcs.crpix
#         final_x, final_y = ax_pix-x_pix, ay_pix+y_pix
        
#         # jpl data
#         x_position, y_position = jpl_horizon_data("17 03 18.38", "-22 05 20.3")
#         ra_hst, dec_hst = x_position - x_pix, y_position + y_pix

#         # Highlight the center pixel of the planet
#         plt.scatter(final_x, final_y, color='red', s=30, label='Mast Observation')
#         plt.scatter(ra_hst, dec_hst, color='blue', s=30, label='JPL observation')
#         draw_circle(ra_hst, dec_hst, 600)
#         plt.title('Mapping Jupiter as seen from HST')
#         plt.xlabel('RA/ Degrees')
#         plt.ylabel('Dec/ Degrees')
#         plt.legend()
#         plt.show()

#     return header

# def time_and_date(header):
#     date_obs = header['DATE-OBS']
#     time_obs = header['TIME-OBS']
#     return date_obs, time_obs 

# def jpl_horizon_data(ra, dec):
#     ra_angle = Angle(ra, unit='hourangle')
#     ra_degrees = ra_angle.deg
#     dec_angle = Angle(dec, unit='deg')
#     dec_degrees = dec_angle.deg
#     return ra_degrees, dec_degrees

# def draw_circle(x, y, radius):
#     circle = plt.Circle((x, y), radius, fill=False, color='yellow', label='Jupiter Border')
#     plt.gca().add_patch(circle)

# file_name = 'j9rlb1i3q_raw.fits'
# header = plot_file(file_name)
# observation_time_date = time_and_date(header)
# print(observation_time_date)

# Figure 2

# from astropy.io import fits
# import matplotlib.pyplot as plt
# from astropy.wcs import WCS
# # from astropy.coordinates import SkyCoord
# from astropy.coordinates import Angle

# def plot_file(file_name):
#     with fits.open(file_name) as hdu_list:
#         header = hdu_list[0].header + hdu_list[1].header
#         w = WCS(header)
#         print(w)
#         fig, ax = plt.subplots()
#         plt.imshow(hdu_list[1].data, origin='lower')

#         # Calculate pixel coordinates based on WCS information
#         ax_pix, ay_pix = w.wcs.crval
#         x_pix, y_pix = w.wcs.crpix
#         final_x, final_y = ax_pix, ay_pix
        
#         # jpl data
#         x_position, y_position = jpl_horizon_data("11 19 13.13","05 39 43.9")
#         ra_hst, dec_hst = x_position , y_position 

#         # Highlight the center pixel of the planet
#         plt.scatter(final_x, final_y, color='red', s=30, label='Mast Observation')
#         plt.scatter(ra_hst, dec_hst, color='blue', s=30, label='JPL observation')
#         draw_circle(ra_hst, dec_hst, 1200)
#         plt.title('Mapping Jupiter as seen from HST')
#         plt.xlabel('RA/ Degrees')
#         plt.ylabel('Dec/ Degrees')
#         plt.legend(loc = "lower right")
#         plt.show()

#     return header

# def time_and_date(header):
#     date_obs = header['DATE-OBS']
#     time_obs = header['TIME-OBS']
#     return date_obs, time_obs 

# def jpl_horizon_data(ra, dec):
#     ra_angle = Angle(ra, unit='hourangle')
#     ra_degrees = ra_angle.deg
#     dec_angle = Angle(dec, unit='deg')
#     dec_degrees = dec_angle.deg
#     return ra_degrees, dec_degrees

# def draw_circle(x, y, radius):
#     circle = plt.Circle((x, y), radius, fill=False, color='yellow', label='Jupiter Border')
#     plt.gca().add_patch(circle)

# file_name = 'ocx841e8q_raw.fits'
# header = plot_file(file_name)
# observation_time_date = time_and_date(header)
# print(observation_time_date)

# Figure 3

from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
# from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle

def plot_file(file_name):
    with fits.open(file_name) as hdu_list:
        header = hdu_list[0].header + hdu_list[1].header
        w = WCS(header)
        print(w)
        fig, ax = plt.subplots()
        plt.imshow(hdu_list[1].data, origin='lower')

        # Calculate pixel coordinates based on WCS information
        ax_pix, ay_pix = w.wcs.crval
        x_pix, y_pix = w.wcs.crpix
        final_x, final_y = (2*x_pix) - ax_pix, (2*y_pix) - ay_pix 
        
        # jpl data
        x_position, y_position = jpl_horizon_data("13 02 46.43" ,"-05 24 24.1")
        ra_hst, dec_hst = (2*x_pix) - x_position, (2*y_pix) - y_position 

        # Highlight the center pixel of the planet
        plt.scatter(final_x, final_y, color='red', s=50, label='Mast Observation')
        plt.scatter(ra_hst, dec_hst, color='blue', s=30, label='JPL observation')
        draw_circle(ra_hst, dec_hst, 1100)
        plt.title('Mapping Jupiter as seen from HST')
        plt.xlabel('RA/ Degrees')
        plt.ylabel('Dec/ Degrees')
        plt.legend(loc = "lower right")
        plt.show()

    return header

def time_and_date(header):
    date_obs = header['DATE-OBS']
    time_obs = header['TIME-OBS']
    return date_obs, time_obs 

def jpl_horizon_data(ra, dec):
    ra_angle = Angle(ra, unit='hourangle')
    ra_degrees = ra_angle.deg
    dec_angle = Angle(dec, unit='deg')
    dec_degrees = dec_angle.deg
    return ra_degrees, dec_degrees

def draw_circle(x, y, radius):
    circle = plt.Circle((x, y), radius, fill=False, color='yellow', label='Jupiter Border')
    return plt.gca().add_patch(circle)

def circle_radius(x_pix, y_pix):
    radius = 


file_name = 'od8k01r0q_raw.fits'
header = plot_file(file_name)
observation_time_date = time_and_date(header)
print(observation_time_date)

# Figure 4

# from astropy.io import fits
# import matplotlib.pyplot as plt
# from astropy.wcs import WCS
# # from astropy.coordinates import SkyCoord
# from astropy.coordinates import Angle

# def plot_file(file_name):
#     with fits.open(file_name) as hdu_list:
#         header = hdu_list[0].header + hdu_list[1].header
#         w = WCS(header)
#         print(w)
#         fig, ax = plt.subplots()
#         plt.imshow(hdu_list[1].data, origin='lower')

#         # Calculate pixel coordinates based on WCS information
#         ax_pix, ay_pix = w.wcs.crval
#         x_pix, y_pix = w.wcs.crpix
#         final_x, final_y = (2*x_pix) - ax_pix, (2*y_pix) - ay_pix 
        
#         # jpl data
#         x_position, y_position = jpl_horizon_data("17 23 09.10"," -22 33 04.7")
#         ra_hst, dec_hst =  2 * x_pix - x_position,  2 * y_pix  - y_position
#         print("center pixel = ", x_pix, y_pix, "jupiter center = ", ax_pix, ay_pix, "final_pix = ", ra_hst, dec_hst)

#         # Highlight the center pixel of the planet
#         plt.scatter(final_x, final_y, color='red', s=30, label='Mast Observation')
#         plt.scatter(ra_hst, dec_hst, color='blue', s=30, label='JPL observation')
#         plt.scatter(ax_pix, ay_pix, color='white', s=20, label='Center Pixel')
#         plt.scatter(x_pix, y_pix, color='green', s=20, label='CRPIX')
#         draw_circle(ra_hst, dec_hst, 1300)
#         plt.title('Mapping Jupiter as seen from HST')
#         plt.xlabel('RA/ Degrees')
#         plt.ylabel('Dec/ Degrees')
#         plt.legend(loc = "lower right")
#         plt.show()

#     return header

# def time_and_date(header):
#     date_obs = header['DATE-OBS']
#     time_obs = header['TIME-OBS']
#     return date_obs, time_obs 

# def jpl_horizon_data(ra, dec):
#     ra_angle = Angle(ra, unit='hourangle')
#     ra_degrees = ra_angle.deg
#     dec_angle = Angle(dec, unit='deg')
#     dec_degrees = dec_angle.deg
#     return ra_degrees, dec_degrees

# def draw_circle(x, y, radius):
#     circle = plt.Circle((x, y), radius, fill=False, color='yellow', label='Jupiter Border')
#     plt.gca().add_patch(circle)

# file_name = 'odxc23qcq_raw.fits'
# header = plot_file(file_name)
# observation_time_date = time_and_date(header)
# print(observation_time_date)


