import os
import numpy as np
import scipy.interpolate as spi
import astropy.units as u


class ThermalEvolution():
    """Class to map thermal evolution parameters from input to output (and vice versa)"""
    def __init__(self, data_directory, data_filename_root, data_filename_suffix, data_legend):
        self.data_directory = data_directory
        self.data_filename_root = data_filename_root
        self.data_filename_suffix = data_filename_suffix
        self.data_legend = data_legend
        self._load_legend()
        self._load_data()

    def _load_legend(self):
        """Load data legend"""
        legend = np.genfromtxt(self.data_legend, dtype='str')
        self.input_parameters = np.unique(np.char.strip(legend, 'T').astype(float), axis=0)
        self.n_simulations = self.input_parameters.shape[0]

    def _load_data(self):
        """Load data"""
        self.thermal_parameters = []
        for simulation_index in self.input_parameters[:, 0].astype(int):
            simulation_string = str(simulation_index).zfill(3)
            data_filepath = os.path.join(self.data_directory, self.data_filename_root + simulation_string + self.data_filename_suffix)
            if os.path.isfile(data_filepath):
                self.thermal_parameters.append(np.genfromtxt(data_filepath))
            else:
                self.thermal_parameters.append(np.ones_like(self.thermal_parameters[-1]) * np.nan)
            #print(self.thermal_parameters[-1].shape)
        assert len(self.thermal_parameters) == self.n_simulations

    def train_interpolator(self, pivot_redshift):
        """Train the interpolator"""
        training_input_parameters = self.input_parameters[:, 1:]

        training_thermal_parameters = []
        training_data_mask = np.ones(self.n_simulations, dtype=np.bool)
        for simulation_number in range(self.n_simulations):
            redshift_index = np.where(self.thermal_parameters[simulation_number][:, 0] == pivot_redshift)[0]
            if redshift_index.shape[0] > 0:
                training_thermal_parameters.append(self.thermal_parameters[simulation_number][redshift_index[0], 1:])
            else:
                training_data_mask[simulation_number] = False
        training_thermal_parameters = np.stack(training_thermal_parameters, axis=0)
        training_input_parameters = training_input_parameters[training_data_mask]
        assert training_input_parameters.shape[0] == training_thermal_parameters.shape[0]

        self._predict_A = spi.Rbf(training_thermal_parameters[:, 0], training_thermal_parameters[:, 1], training_thermal_parameters[:, 2], training_input_parameters[:, 0])
        self._predict_B = spi.Rbf(training_thermal_parameters[:, 0], training_thermal_parameters[:, 1],
                                 training_thermal_parameters[:, 2], training_input_parameters[:, 1])

    def predict_A(self, T0, gamma, filtering_length):
        """Predict the value of A (multiplicative correction to photoheating rates) required to give the given thermal
        parameters (T0 in K, gamma and filtering_length in ckpc) at the chosen pivot redshift"""
        return self._predict_A(T0.to(u.K), gamma, filtering_length.to(u.kpc))

    def predict_B(self, T0, gamma, filtering_length):
        """Predict the value of B (exponent of density-dependent correction to photoheating rates) required to give the
        given thermal parameters (T0 in K, gamma and filtering_length in ckpc) at the chosen pivot redshift"""
        return self._predict_B(T0.to(u.K), gamma, filtering_length.to(u.kpc))
