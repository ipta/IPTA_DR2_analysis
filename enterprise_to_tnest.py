#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 17:04:04 2019

@author: dreardon

enterprise_to_tnest.py: converts enterprise json to temponest parameters
"""

import json
import numpy as np
import glob


def write_params(psrname, line):
    with open(pars_dir + '/' + psrname + '.par', "a+") as myfile:
        myfile.write(line + '\n')

# Enter directory names for parameters and noise files
pars_dir = '/your/par/files/directory'
noise_dir = '/your/enterprise/noise/files/directory'

for noisefile in sorted(glob.glob(noise_dir + '*.json')):
    
    with open(noisefile) as json_file:
        params = json.load(json_file)

    param_count = 0
    dmexp_count = 0
    magnexp_count = 0
    for key, value in params.items():
        if key == 'nmodel':
            continue
        psrname = key.split('_')[0]
        if param_count == 0:
            write_params(psrname, ' ')
            write_params(psrname, '### Noise parameters measured with enterprise')
        param_count += 1
        key = key.replace(psrname+'_', '')
        print(param_count, ' ', key)
        if '_efac' in key:
            key = key.replace('_efac', '')
            line = 'TNEF -group ' + key + ' ' + str(value)
            write_params(psrname, line)
        elif '_log10_equad' in key:
            key = key.replace('_log10_equad', '')
            line = 'TNEQ -group ' + key + ' ' + str(value)
            write_params(psrname, line)
        elif 'dmexp' in key:
            if 'log10_Amp' in key:
                dmexp_count += 1 
                write_params(psrname, 'EXPINDEX_{0} -2'.format(int(dmexp_count)))
                write_params(psrname, 'EXPPH_{0} {1}'.format(int(dmexp_count), 10**value))
            elif 'log10_tau' in key:
                write_params(psrname, 'EXPTAU_{0} {1}'.format(int(dmexp_count), 10**value))
            elif '_t0' in key:
                write_params(psrname, 'EXPEP_{0} {1}'.format(int(dmexp_count), value))
        elif 'magnexp' in    key:
               if 'idx' in key:
                   magnexp_count += 1
                   write_params(psrname, 'EXPINDEX_{0} {1}'.format(int(magnexp_count), -value))
               elif 'log10_Amp' in key:
                   write_params(psrname, 'EXPPH_{0} {1}'.format(int(magnexp_count), 10**value))
               elif 'log10_tau' in    key:
                      write_params(psrname, 'EXPTAU_{0} {1}'.format(int(magnexp_count), 10**value))
               elif '_t0' in key:
                      write_params(psrname, 'EXPEP_{0} {1}'.format(int(magnexp_count), value))
        elif 'dm_gp' in key:
            if 'gamma' in key:
                line = 'TNDMGam ' + str(value)
                write_params(psrname, line)
            elif 'log10_A' in key:
                line = 'TNDMAmp ' + str(np.log10(np.sqrt(12*np.pi**2)*10**value))
                write_params(psrname, line)
                line = 'TNDMC 90'
                write_params(psrname, line)
                line = 'TNSubtractDM 1'
                write_params(psrname, line)
        elif 'red_noise' in key:
            if 'low_freq' in key:
                if 'gamma' in key:
                    linelow = 'TNBandNoise 0 1000 ' + 'AMP ' + str(np.log10(np.sqrt(12*np.pi**2)*10**value)) + ' 90'
                elif 'log10_A' in key:
                    write_params(psrname, linelow.replace('AMP', str(value)))
            elif 'mid_freq' in key:
                if 'gamma' in key:
                    linemid = 'TNBandNoise 1000 2000 ' + 'AMP ' + str(np.log10(np.sqrt(12*np.pi**2)*10**value)) + ' 90'
                elif 'log10_A' in key:
                    write_params(psrname, linemid.replace('AMP', str(value)))
                    write_params(psrname, 'TNRedGam 2')
                    write_params(psrname, 'TNRedAmp -20')
                    write_params(psrname, 'TNRedC 2')
                    write_params(psrname, 'TNSubtractRed 1')
            elif 'high_freq' in key:
                if 'gamma' in key:
                    linehigh = 'TNBandNoise 2000 5000 ' + 'AMP ' + str(np.log10(np.sqrt(12*np.pi**2)*10**value)) + ' 90'
                elif 'log10_A' in key:
                    write_params(psrname, linehigh.replace('AMP', str(value)))
            else:
                if 'gamma' in key:
                    line = 'TNRedGam ' + str(value)
                    write_params(psrname, line)
                elif 'log10_A' in key:
                    line = 'TNRedAmp ' + str(value)
                    write_params(psrname, line)
                    line = 'TNRedC 90'
                    write_params(psrname, line)
                    line = 'TNSubtractRed 1'
                    write_params(psrname, line)
    print('Wrote temponest parameters for {0}'.format(psrname))
