"""Generates a Latin hypercube, and generates corner plot."""

import sys
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import configargparse
# our modules below
import camb_cosmo
import fit_linP
import sim_params_cosmo
import sim_params_space
import read_genic
import write_config
import latin_hypercube

# get options from command line
parser = configargparse.ArgumentParser()
parser.add_argument('-c', '--config', required=False, is_config_file=True, help='config file path')
parser.add_argument('--basedir', type=str, help='Base directory where all sims will be stored (crashes if it already exists)',required=True)
parser.add_argument('--add_slope', action='store_true', help='Add parameter describing slope of linear power',required=False)
parser.add_argument('--add_running', action='store_true', help='Add parameter describing running of linear power',required=False)
parser.add_argument('--add_mu_H', action='store_true', help='Add parameter to boost heating in Hydrogen',required=False)
parser.add_argument('--nsamples', type=int, default=10, help='Number of samples in Latin hypercube')
parser.add_argument('--ngrid', type=int, default=64, help='Number of particles per side in simulation')
parser.add_argument('--box_Mpc', type=float, default=50.0, help='Simulation box size (in Mpc; no h normalisation)')
parser.add_argument('--zs', type=str, help='Comma-separated list of redshifts (including last snapshot)')
parser.add_argument('--z_star', type=float, default=3.0, help='Pivot redshift')
parser.add_argument('--kp_Mpc', type=float, default=0.7, help='Pivot wavenumber (in Mpc; no h normalisation)')
parser.add_argument('--seed', type=int, default=123, help='Random seed to setup Latin hypercube')
parser.add_argument('--verbose', action='store_true', help='Print runtime information',required=False)
args = parser.parse_args()

print('--- print options from parser ---')
print(args)
print("----------")
print(parser.format_help())
print("----------")
print(parser.format_values()) 
print("----------")

verbose=args.verbose

# transform input string to list of sorted redshifts
if args.zs:
    zs=np.sort([float(z) for z in args.zs.split(',')])[::-1]
    print('will use input redshifts',zs)
else:
    zs=None

# setup parameter space
param_space=sim_params_space.SimulationParameterSpace(filename=args.config,
                    add_slope=args.add_slope,add_running=args.add_running,
                    add_mu_H=args.add_mu_H, z_star=args.z_star, kp_Mpc=args.kp_Mpc)
params=param_space.params

# get pivot point
z_star=params['Om_star']['z_star']
kp_Mpc=params['Delta2_star']['kp_Mpc']

# print parameter information
if verbose:
    print('z_star =',z_star)
    print('kp_Mpc =',kp_Mpc)
    for key,param in params.items():
        print(key,param)

# get parameter ranges
Npar=len(params)
param_limits=np.empty([Npar,2])
for key,param in params.items():
    ip=param['ip']
    param_limits[ip][0]=param['min_val']
    param_limits[ip][1]=param['max_val']

# generate Latin hypercube 
nsamples=args.nsamples
seed=args.seed
cube=latin_hypercube.get_hypercube_samples(param_limits, nsamples, 
        prior_points = None, seed=seed)

# print information about cube
if verbose:
    print('# samples =',nsamples)
    print('random seed',seed)
    print('initial points in cube')
    print(cube)

# get fiducial cosmology
cosmo_fid = camb_cosmo.get_cosmology()
if verbose:
    camb_cosmo.print_info(cosmo_fid)

# setup fiducial linear power model
linP_model_fid=fit_linP.LinearPowerModel(cosmo=cosmo_fid,z_star=z_star,
            k_units='Mpc',kp=kp_Mpc)
if verbose:
    print('fiducial linear power parameters',linP_model_fid.get_params())

# make sure the base directory does not exist
basedir=args.basedir
if os.path.exists(basedir):
    raise ValueError(basedir+' already exists')
os.mkdir(basedir)

# write file with description of the hypercube
write_config.write_cube_json_file(basedir,params,cube)
for sample in range(nsamples):
    sim_params=cube[sample]
    if verbose:
        print(sample,sim_params)
    cosmo_sim=sim_params_cosmo.cosmo_from_sim_params(params,sim_params,
            linP_model_fid,verbose=verbose)

    #nCDM parameters
    if 'Alpha' in params:
        ip_alpha = params['Alpha']['ip']
        alpha = sim_params[ip_alpha]
    else:
        alpha = 0.

    if 'Beta' in params:
        ip_beta = params['Beta']['ip']
        beta = sim_params[ip_beta]
    else:
        beta = 1.

    if 'Gamma' in params:
        ip_gamma = params['Gamma']['ip']
        gamma= sim_params[ip_gamma]
    else:
        gamma = 0.

    # figure out heating amplitude for Helium and Hydrogen
    if 'mu_He' in params:
        ip_He=params['mu_He']['ip']
        mu_He=sim_params[ip_He]
    else:
        mu_He=1.0
    if 'mu_H' in params:
        ip_H=params['mu_H']['ip']
        mu_H=sim_params[ip_H]
    else:
        mu_H=1.0

    sim_dir=basedir+'/sim_'+str(sample)+'/'
    os.mkdir(sim_dir)

    # write GenIC and MP-Gadget parameters, for both simulations in pair
    if verbose:
        print('write config files for GenIC and Gadget')
    write_config.write_genic_file(sim_dir,cosmo_sim, alpha=alpha, beta=beta, gamma=gamma,
            Ngrid=args.ngrid,box_Mpc=args.box_Mpc,paired=False)
    zs=write_config.write_gadget_file(sim_dir,cosmo_sim,mu_H=mu_H,mu_He=mu_He,
            Ngrid=args.ngrid,zs=zs)

    # construct linear power model and store in JSON format
    linP_model_sim=fit_linP.LinearPowerModel(cosmo=cosmo_sim,z_star=z_star,
            k_units='Mpc',kp=kp_Mpc)

    if verbose:
        print('write JSON file for simulation pair')
    write_config.write_sim_json_file(sim_dir,params,sim_params,linP_model_sim,
            zs=zs)

print('finished')
