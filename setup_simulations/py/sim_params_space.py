"""Dictionary describing the simulation parameter space."""

import numpy as np

class SimulationParameterSpace(object):
    """Describe simulation parameter space, to be used by latin hypercube."""

    def __init__(self,filename=None,add_slope=True,add_running=False,
                add_mu_He=True,add_mu_H=False, z_star=4., kp_Mpc=2.):
        """Construct space from file, or using default setting"""

        if filename is None:
            self._default_setup_nCDM(add_slope,add_running,add_mu_He,add_mu_H, z_star=z_star, kp_Mpc=kp_Mpc)
        else:
            self._setup_from_file(filename,add_slope,add_running,
                            add_mu_He,add_mu_H, z_star=z_star, kp_Mpc=kp_Mpc)


    def _default_setup(self,add_slope,add_running,add_mu_He,add_mu_H, z_star=3., kp_Mpc=0.7):
        """Default setup of parameter space"""

        params={}
        params['Om_star']={'ip':len(params), 'min_val':0.955, 'max_val':0.975,
                'z_star':z_star, 'latex':r'$\Omega_\star$'}
        params['Delta2_star']={'ip':len(params), 'min_val':0.25, 'max_val':0.45,
                'z_star':z_star, 'kp_Mpc':kp_Mpc, 'latex':r'$\Delta^2_\star$'}
        if add_slope:
            params['n_star']={'ip':len(params), 'min_val':-2.35,
                    'max_val':-2.25, 'z_star':z_star, 'kp_Mpc':kp_Mpc,
                    'latex':r'$n_\star$'}
        if add_running:
            params['alpha_star']={'ip':len(params), 'min_val':-0.265,
                    'max_val':-0.165, 'z_star':z_star, 'kp_Mpc':kp_Mpc,
                    'latex':r'$\alpha_\star$'}
        if add_mu_He:
            params['mu_He']={'ip':len(params), 'min_val':0.5, 'max_val':2.0,
                    'latex':r'$\mu_{\rm He}$'}
        if add_mu_H:
            params['mu_H']={'ip':len(params), 'min_val':0.5, 'max_val':2.0,
                    'latex':r'$\mu_{\rm H}$'}
        self.params=params


    def _default_setup_nCDM(self,add_slope,add_running,add_mu_He,add_mu_H, z_star=4., kp_Mpc=2.):
        """Default setup of parameter space"""

        params={}
        params['Om_star']={'ip':len(params), 'min_val':0.977, 'max_val':0.986,
                'z_star':z_star, 'latex':r'$\Omega_\star$'}
        params['Delta2_star']={'ip':len(params), 'min_val':0.35, 'max_val':0.55,
                'z_star':z_star, 'kp_Mpc':kp_Mpc, 'latex':r'$\Delta^2_\star$'}
        params['Alpha'] = {'ip': len(params), 'min_val': 0., 'max_val': 0.1, 'latex': r'$\alpha$'}
        params['Beta'] = {'ip': len(params), 'min_val': 0., 'max_val': 10., 'latex': r'$\beta$'}
        params['Gamma'] = {'ip': len(params), 'min_val': -10., 'max_val': 0., 'latex': r'$\gamma$'}
        if add_slope:
            params['n_star']={'ip':len(params), 'min_val':-2.53,
                    'max_val':-2.43, 'z_star':z_star, 'kp_Mpc':kp_Mpc,
                    'latex':r'$n_\star$'}
        if add_running:
            params['alpha_star']={'ip':len(params), 'min_val':-0.265,
                    'max_val':-0.165, 'z_star':z_star, 'kp_Mpc':kp_Mpc,
                    'latex':r'$\alpha_\star$'}
            print('Prior limits on alpha_star could be reconsidered!')
        if add_mu_He:
            params['mu_He']={'ip':len(params), 'min_val':0.5, 'max_val':2.0,
                    'latex':r'$\mu_{\rm He}$'}
        if add_mu_H:
            params['mu_H']={'ip':len(params), 'min_val':0.5, 'max_val':2.0,
                    'latex':r'$\mu_{\rm H}$'}
        self.params=params


    def _setup_from_file(self,filename,add_slope,add_running,
                    add_mu_He,add_mu_H, z_star=3., kp_Mpc=0.7):
        print('should implement setup from file')
        self._default_setup(add_slope,add_running,add_mu_He,add_mu_H, z_star=z_star, kp_Mpc=kp_Mpc)
        #raise ValueError('implement setup_from_file')
