import os
import numpy as np
import camb
import camb_cosmo
import fit_linP

class ReconstructedCosmology(object):
    """Given fiducial cosmology, and set of linear power parameters, 
        reconstruct a cosmology object."""

    def __init__(self,zs,linP_model=None,cosmo_fid=None,use_constant_f=False,
                    verbose=False):
        """Setup from linear power model and fiducial cosmology, as well
            as redshifts that we will want to evaluate (zs)."""

        self.verbose=verbose
        # store redshifts that will be evaluated
        self.zs=zs

        # fiducial cosmology
        if cosmo_fid:
            if self.verbose: print('use input fiducial cosmology')
            self.cosmo_fid=cosmo_fid
        else:
            if self.verbose: print('use default fiducial cosmology')
            self.cosmo_fid=camb_cosmo.get_cosmology()

        # compute CAMB results for fiducial cosmology
        self.results_fid=camb.get_results(self.cosmo_fid)

        # input model describing linear power around z_star and kp_kms
        if linP_model:
            if self.verbose: print('use input linear power model')
            self.linP_model=linP_model
            assert linP_model.k_units is 'kms', 'input linP_model not in kms'
            # compute linear power model for fiducial cosmology
            self.z_star=linP_model.z_star
            self.kp_kms=linP_model.kp
            self.linP_model_fid=fit_linP.LinearPowerModel(cosmo=self.cosmo_fid,
                            z_star=self.z_star,k_units='kms',kp=self.kp_kms)
        else:
            if self.verbose: print('no input linear power model')
            self.linP_model=None
            # compute linear power model for fiducial cosmology
            self.linP_model_fid=fit_linP.LinearPowerModel(cosmo=self.cosmo_fid)
            self.z_star=self.linP_model_fid.z_star
            self.kp_kms=self.linP_model_fid.kp


        # whether to model z-evolution of logarithmic growth rate with fiducial
        self.use_constant_f=use_constant_f

        # get Hubble at z_star for fiducial cosmology, used to compute kp_Mpc
        self.H_star_fid=self.results_fid.hubble_parameter(self.z_star)
        self.dkms_dMpc_star_fid=self.H_star_fid/(1+self.z_star)
        self.kp_Mpc=self.kp_kms*self.dkms_dMpc_star_fid

        # store Hubble parameter at all redshifts for fiducial cosmology
        self._compute_Hz_fid()

        # store linear power at all redshifts for fiducial cosmology
        self._compute_f_p_fid()

        # store linear power at all redshifts for fiducial cosmology
        self._compute_linP_Mpc_fid()


    def _compute_Hz_fid(self):
        """ Compute Hubble parameter in fiducial cosmology, at all redshifts."""

        self.Hz_fid=[]
        for z in self.zs:
            Hz = self.results_fid.hubble_parameter(z)
            self.Hz_fid.append(Hz)


    def _compute_f_p_fid(self):
        """ Compute growth rate in fiducial cosmology, at all redshifts. """

        self.f_p_fid=[]
        for z in self.zs:
            f_p = fit_linP.compute_f_star(self.cosmo_fid,z_star=z,
                                                kp_Mpc=self.kp_Mpc)
            self.f_p_fid.append(f_p)


    def _compute_linP_Mpc_fid(self):
        """ Compute linear power in fiducial cosmology, at all redshifts. """

        # call CAMB to get linear power for fiducial cosmology, in Mpc
        k_Mpc, zs_out, P_Mpc = camb_cosmo.get_linP_Mpc(self.cosmo_fid,self.zs)
        # make sure we didn't change the order of the redshift outputs
        assert zs_out[0] == self.zs[0], 'CAMB redshifts not sorted'

        self.k_Mpc = k_Mpc
        self.linP_Mpc_fid = P_Mpc


    def get_linP_Mpc_params(self,linP_model=None):
        """Reconstruct linear power (in Mpc) for input linP_model and fit
            linear power parameters at each redshift. """

        if linP_model:
            if self.verbose: print('use input linP_model')
            self.linP_model=linP_model
        else:
            if self.verbose: print('use stored linP_model')
            if not self.linP_model:
                raise ValueError('we need a model in get_linP_Mpc_Params')

        # wavenumbers that will be used in fit
        kp_Mpc=self.kp_Mpc
        kmin_Mpc=0.5*kp_Mpc
        kmax_Mpc=2.0*kp_Mpc
        xmin=kmin_Mpc/kp_Mpc
        xmax=kmax_Mpc/kp_Mpc
        x=self.k_Mpc/kp_Mpc

        linP_Mpc_params=[]
        for iz,z in enumerate(self.zs):
            # get information from fiducial cosmology at this redshift
            linP_Mpc_fid=self.linP_Mpc_fid[iz]
            # reconstruct logarithmic growth rate at the redshift
            f_p=self.reconstruct_f_p_iz(iz)
            # reconstruct linear power at the redshift (in Mpc)
            linP_Mpc=self.reconstruct_linP_Mpc(z,linP_Mpc_fid=linP_Mpc_fid)
            # fit polynomial describing log linear power
            linP_fit=fit_linP.fit_polynomial(xmin,xmax,x,linP_Mpc,deg=2)
            # compute parameters used in emulator
            lnA_p=linP_fit[0]
            Delta2_p=np.exp(lnA_p)*kp_Mpc**3/(2*np.pi**2)
            n_p=linP_fit[1]
            # note that the curvature is alpha/2
            alpha_p=2.0*linP_fit[2]
            params={'Delta2_p':Delta2_p,'n_p':n_p,'alpha_p':alpha_p,'f_p':f_p}
            linP_Mpc_params.append(params)

        return linP_Mpc_params


    def reconstruct_linP_Mpc(self,z,linP_Mpc_fid):
        """ Use fiducial cosmology and linP_params to reconstruct power (Mpc)"""

        # pivot points
        z_star=self.z_star
        kp_kms=self.kp_kms

        # reconstruct starting from the power in the fiducial cosmology
        linP_Mpc = np.copy(linP_Mpc_fid)

        # get parameters describing linear power for input cosmology
        f_star=self.linP_model.get_f_star()
        linP_kms_params=self.linP_model.linP_params['linP_kms']
        # get parameters describing linear power for fiducial cosmology
        f_star_fid=self.linP_model_fid.get_f_star()
        linP_kms_params_fid=self.linP_model_fid.linP_params['linP_kms']
        # get relative parameters
        df_star=f_star-f_star_fid

        # modify shape based on fits done in velocity units
        # B(q) = P(z_star,q) / P_fid(z_star,q)
        # linP_kms_params actually store log(B) vs log(k_kms/kp_kms)
        lnk_kp = np.log(self.k_Mpc / self.kp_Mpc)
        lnB = linP_kms_params(lnk_kp)-linP_kms_params_fid(lnk_kp)
        linP_Mpc *= np.exp(lnB)

        # correct for linear growth with respect to fiducial cosmology
        # [D(z) / D(z_star)] / [D_fid(z) / D_fid(z_star)]
        D_correct = 1-df_star*(z-z_star)/(1+z_star)
        linP_Mpc *= D_correct**2

        return linP_Mpc


    def reconstruct_Hubble_iz(self,iz):
        """ Use fiducial cosmology and g_star to reconstruct Hubble parameter"""

        Hz_fid=self.Hz_fid[iz]
        z=self.zs[iz]
        return self.reconstruct_Hubble(z,Hz_fid=Hz_fid)


    def reconstruct_Hubble(self,z,Hz_fid=None):
        """ Use fiducial cosmology and g_star to reconstruct Hubble parameter"""

        if not Hz_fid:
            Hz_fid=self.results_fid.hubble_parameter(z)
        # compute difference in acceleration
        g_star=self.linP_model.get_g_star()
        g_star_fid=self.linP_model_fid.get_g_star()
        # compute Hubble parameter in input cosmology
        z_star=self.z_star
        Hz = Hz_fid * (1+3/2*(g_star-g_star_fid)*(z-z_star)/(1+z_star))
        return Hz


    def reconstruct_f_p_iz(self,iz):
        """ Use fiducial cosmology and f_star to reconstruct logarithmic
            growth rate f (around kp_Mpc)"""

        f_p_fid=self.f_p_fid[iz]
        z=self.zs[iz]
        return self.reconstruct_f_p(z,f_p_fid=f_p_fid)


    def reconstruct_f_p(self,z,f_p_fid=None):
        """ Use fiducial cosmology and f_star to reconstruct logarithmic
            growth rate f (around kp_Mpc)"""

        if self.use_constant_f:
            return self.linP_model.get_f_star()

        # compute f in fiducial cosmology
        if not f_p_fid:
            f_p_fid=fit_linP.compute_f_star(self.cosmo_fid,z_star=z,
                            kp_Mpc=self.kp_Mpc)
        # correct using difference in f_star
        f_star=self.linP_model.get_f_star()
        f_star_fid=self.linP_model_fid.get_f_star()
        df_star=f_star-f_star_fid
        f_p=f_p_fid+(f_star-f_star_fid)
        return f_p

