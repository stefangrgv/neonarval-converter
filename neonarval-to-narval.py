# NeoNarval -> Narval converter
#
# Converts NeoNarval observations (.fits files) into the format of its predecessor, Narval (.s)
#
# Stefan Georgiev
# stefan.grgv@gmail.com
#
# 07 July 2020

from tkinter.filedialog import askopenfilenames
import tkinter as tk
import numpy as np
from astropy.io import fits
import datetime as dt
import os
import sys

##########################################################################################################################
# FUNCTIONS

# doppler shift formula
c = 299792.458
def rvhelio(wl, berv):
    return wl + wl*(berv/c)

# read .fits and return contents 
def get_fits_data(image_file, n, N):
    print('Loading ' + image_file.split('/')[-1] + '... (file %i/'%(n+1) + '%i)'%N)

    a = fits.open(image_file)
    wl=a[1].data['Wavelength1']
    I=a[1].data['Intensity']
    Stk=np.nan_to_num(a[1].data['Polarization'])
    error=np.nan_to_num(a[1].data['Error'])
    null=np.nan_to_num(a[1].data['Null'])
    header=a[0].header
    a.close()

    star    =   header['OBJECT'].replace(' ', '').upper().replace(' ', '').replace('-', '').replace('_', '')

    if desired_stars:
        if star not in desired_stars:
            print('Skipped: object ' + star + ' not in the list of desired stars')
            return None

    stokes  =   header['STOKNAME']

    # I found that for some of the data the DATE-FIT entry in the header is in a different format; the following block allows the script to still convert such observations
    try:
        date_format_ok = True
        datefit =   dt.datetime.strptime(header['DATE-FIT'][:-4], '%Y-%m-%dT%H:%M:%S')
    except ValueError:
        date_format_ok = False
        datefit =   header['DATE-FIT']

    if date_format_ok:
        date_obs = datefit.date()

        if datefit.hour < 12:
            date_obs = date_obs - dt.timedelta(1)

        date_obs_str = date_obs.strftime('%d%b%-y').lower()
    else:
        date_obs_str = datefit


    return(wl, I, Stk, error, null, header, star, stokes, date_obs_str, datefit, date_format_ok)

# do the actual conversion
def convert(wl, I, Stk, error, null, header, star, stokes, date_obs_str, datefit, date_format_ok):
    if date_format_ok and not os.path.exists(date_obs_str):
        os.makedirs(date_obs_str)

    # find out the sequence number of the current observation; note that this is with respect to the current file structure and not to the observations themselves
    if date_format_ok:
        sequence_number = 1
        while True:
            output_file = star.lower() + '_' + stokes.lower() + '_%02i'%(sequence_number)
            output_file = os.path.join(date_obs_str, output_file)
          
            if not os.path.exists(output_file + '.s'):
                break

            sequence_number += 1

    else:
        output_file = star.lower() + '_' + stokes.lower() + '_' + datefit

    print('Writing ' + output_file + '...')

    # locate echelle orders in the spectrum
    start_of_order, end_of_order = [0], []
    for i in range(1, len(wl)):
        if wl[i-1] > wl[i]:
            end_of_order.append(i-1)
            start_of_order.append(i)
    end_of_order.append(len(wl)-1)

    orders = []
    for i in range(len(end_of_order)):
        orders.append([start_of_order[i], end_of_order[i]])

    
    # angstroms to nanometers    
    wl /= 10

    # transform the data into the heliocentric restframe
    berv = header['BERV']
    wl = rvhelio(wl, berv)

    # orders are reversed in neonarval and the lambda decreases from one order to the next: we want to reverse that
    orders = reversed(orders)

    # extract only this wavelength range: LSD breaks if you extract everything
    wl_start, wl_end = 380, 1050
    good_wl = []
    for o in orders:
        start_tmp = o[0]
        end_tmp   = o[1] + 1
        for i in range(start_tmp, end_tmp):
            if wl[i] > wl_start and wl[i] < wl_end:
                good_wl.append(i) # we're making a list of good lambda indices to use in the next block


    # write the .s    
    wl_lines = len(good_wl)
    with open(output_file + '.s', 'w') as f:
        f.write('***Reduced spectrum of \'{0}\'\n'.format(star))
        f.write(' {0} {1}\n'.format(wl_lines,5))
        for i in good_wl:
            f.write('{:10.4f} '.format(wl[i]))
            f.write('{: 8.4e} '.format(I[i]))
            f.write('{: 8.4e} '.format(Stk[i]))
            f.write('{: 8.4e} '.format(null[i]))
            f.write('{: 8.4e} '.format(null[i]))
            f.write('{: 8.4e}\n'.format(error[i]))

    # write the .header
    with open(output_file + '.header', 'w') as f:
        for keyword, data in header.items():
            f.write(str(keyword) + ':\t' + str(data) + '\n')

##########################################################################################################################
# MAIN PROGRAM

# handle user-selected stars, if any
desired_stars = []
for j in range(len(sys.argv)):
    if j > 0:
        desired_stars.append(sys.argv[j].replace('-', '').replace('_', '').replace(',', '').upper())

# select one or more neonarval .fits
tk.Tk().withdraw() # get rid of an extra window
neonarval_input = askopenfilenames(title='Select NeoNarval observations', filetypes=[('NeoNarval fits', '*.fits *.FITS')])
if not neonarval_input:
    print('No files selected.')
    sys.exit()

# convert the selected observations into .s
N = len(neonarval_input)
for n, image_file in enumerate(neonarval_input):
    fits_data = get_fits_data(image_file, n, N)
    if fits_data: # fits_data will be None if the current object is not in the desired_stars list
        wl, I, Stk, error, null, header, star, stokes, date_obs_str, datefit, date_format_ok = fits_data
        convert(wl, I, Stk, error, null, header, star, stokes, date_obs_str, datefit, date_format_ok)

