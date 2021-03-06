import numpy as np
import sys
import os
import json
import scipy.interpolate
import p1d_arxiv
import poly_p1d


class LinearEmulator(object):
    """Linear interpolation emulator for flux P1D."""

    def __init__(self,basedir='../mini_sim_suite/',
            p1d_label='p1d',skewers_label='Ns50_wM0.1',
            emulate_running=False,emulate_growth=False,
            deg=4,kmax_Mpc=10.0,max_arxiv_size=None,verbose=True):
        """Setup emulator from base sim directory and label identifying skewer
            configuration (number, width)"""

        self.verbose=verbose

        # read all files with P1D measured in simulation suite
        self.arxiv=p1d_arxiv.ArxivP1D(basedir,p1d_label,skewers_label,
                    max_arxiv_size,verbose)

        # for each model in arxiv, fit smooth function to P1D
        self._fit_p1d_in_arxiv(deg,kmax_Mpc)

        # setup parameter space to be used in emulator
        self._setup_param_space(emulate_running,emulate_growth)

        # for each order in polynomial, setup interpolation object
        self._setup_interp(deg)
        

    def _fit_p1d_in_arxiv(self,deg,kmax_Mpc):
        """For each entry in arxiv, fit polynomial to log(p1d)"""
        
        for entry in self.arxiv.data:
            k_Mpc = entry['k_Mpc']
            p1d_Mpc = entry['p1d_Mpc']
            fit_p1d = poly_p1d.PolyP1D(k_Mpc,p1d_Mpc,kmin_Mpc=1.e-3,
                    kmax_Mpc=kmax_Mpc,deg=deg)
            entry['fit_p1d'] = fit_p1d


    def _setup_param_space(self,emulate_running,emulate_growth):
        """Set order of parameters in emulator"""

        self.params=['Delta2_p','n_p']
        if emulate_running:
            self.params.append('alpha_p')
        if emulate_growth:
            self.params.append('growht')
        self.params += ['mF','sigT_Mpc','gamma']
        if self.verbose:
            print('parameter names in emulator',self.params)


    def _setup_interp(self,deg):
        """For each order in polynomial, setup interpolation object"""

        # for each parameter in params, get values from arxiv
        point_params=[]
        for par in self.params:
            values = np.array([entry[par] for entry in self.arxiv.data])
            point_params.append(values)
        self.points=np.vstack(point_params).transpose()

        N=len(self.arxiv.data)
        self.linterps=[]
        for p in range(deg+1):
            print('setup interpolator for coefficient',p)
            values = [entry['fit_p1d'].lnP[p] for entry in self.arxiv.data] 
            linterp = scipy.interpolate.LinearNDInterpolator(self.points,values)
            self.linterps.append(linterp)
            # it is good to try the interpolator to finish the setup
            # (it might help to avoid thread issues later on)
            test_point=np.median(self.points,axis=0)
            print(test_point,'test',linterp(test_point))


    def _point_from_model(self,model):
        """Extract model parameters from dictionary in the right order"""

        point=[]
        for par in self.params:
            point.append(model[par])

        return np.array(point)


    def emulate_p1d_Mpc(self,model,k_Mpc):
        """Return emulate 1D power spectrum at input k values"""

        if self.verbose: print('asked to emulate model',model)

        # get interpolation point from input model
        point = self._point_from_model(model)
        if self.verbose: print('evaluate point',point)

        # emulate coefficients for PolyP1D object (note strange order of coeffs)
        Npar=len(self.linterps)
        coeffs=np.empty(Npar)
        for i in range(Npar):
            coeffs[Npar-i-1] = self.linterps[i](point)
        if self.verbose: print('got coefficients',coeffs)

        # set P1D object
        kmin_Mpc=self.arxiv.data[0]['fit_p1d'].kmin_Mpc
        smooth_p1d = poly_p1d.PolyP1D(lnP_fit=coeffs,kmin_Mpc=kmin_Mpc)

        return smooth_p1d.P_Mpc(k_Mpc)

