import numpy as np
import thermal_model


def emulate_p1d(z,emu,dkms_dMpc,mf_model,T_model,linP_Mpc_params):
    """Emulate 1D power given model and redshift"""
    
    # get emulator parameters for linear power, at this redshift (in Mpc)
    model=linP_Mpc_params
    # get emulator parameters for nuisance models, at this redshift
    model['mF']=mf_model.get_mean_flux(z)
    model['gamma']=T_model.get_gamma(z)
    T0=T_model.get_T0(z)
    sigT_kms=thermal_model.thermal_broadening_kms(T0)
    model['sigT_Mpc']=sigT_kms/dkms_dMpc

    return emu.emulate_p1d(model)


def get_chi2(data,cosmo_fid,emu,rec_cosmo,mf_model,T_model,
            linP_Mpc_params=None):
    """Compute chi2 given data, fiducial cosmology, emulator, 
        reconstructed cosmology, nuisance params and linear power params."""

    # while testing, use only a handful of bins
    zs=data.z[::4]
    
    # check if linear power parameters have been cached
    if linP_Mpc_params is None:
        linP_Mpc_params=rec_cosmo.get_linP_Mpc_params(zs)
    
    Nz=len(zs)
    chi2=0
    for iz in range(Nz):
        # acess data for this redshift
        z=zs[iz]
        # get conversion from Mpc to km/s
        dkms_dMpc=rec_cosmo.reconstruct_Hubble(z)/(1+z)
        # get data
        p1d=data.get_Pk_iz(iz)
        cov=data.get_cov_iz(iz)
        emu_k_Mpc, emu_p1d_Mpc = emulate_p1d(z,emu,dkms_dMpc,mf_model,T_model,
                    linP_Mpc_params[iz])
        # translate to km/s
        emu_k_kms = emu_k_Mpc / dkms_dMpc
        emu_P_kms = emu_p1d_Mpc * dkms_dMpc
        # interpolate to wavenumbers in data
        k_kms=data.k
        emu_p1d = np.interp(k_kms,emu_k_kms,emu_P_kms)
        # compute chi2 for this redshift bin
        icov = np.linalg.inv(cov)
        diff = (p1d-emu_p1d)
        chi2_z = np.dot(np.dot(icov,diff),diff)
        chi2 += chi2_z
        
    return chi2