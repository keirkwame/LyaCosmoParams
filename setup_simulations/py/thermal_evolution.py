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

    def train_interpolator(self, pivot_redshift, use_parameter=[True, True, True]):
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

        #training_thermal_parameters[:, 1] *= 1.e+4
        #training_thermal_parameters[:, 2] *= 1.e+2
        training_thermal_parameters = training_thermal_parameters[:, use_parameter]
        self.training_thermal_parameters_pivot_redshift = training_thermal_parameters
        self._minimum_training_thermal_parameters = np.min(training_thermal_parameters, axis=0)
        training_thermal_parameters_zeroed = training_thermal_parameters - self._minimum_training_thermal_parameters[np.newaxis, :]
        self._maximum_training_thermal_parameters_zeroed = np.max(training_thermal_parameters_zeroed, axis=0)
        training_thermal_parameters = training_thermal_parameters_zeroed / self._maximum_training_thermal_parameters_zeroed[np.newaxis, :]
        assert np.min(training_thermal_parameters) == 0.
        assert np.max(training_thermal_parameters) == 1.
        training_thermal_parameters_list = np.split(training_thermal_parameters, training_thermal_parameters.shape[1], axis=1)

        self._predict_A = spi.Rbf(*training_thermal_parameters_list, training_input_parameters[:, 0])
        self._predict_B = spi.Rbf(*training_thermal_parameters_list, training_input_parameters[:, 1])

    def map_thermal_parameters_list_to_unit_hypercube(self, thermal_parameters_list):
        """Map thermal parameters to unit hypercube (limits determined by space spanned by training data)"""
        for i in range(len(thermal_parameters_list)):
            thermal_parameters_list[i] = thermal_parameters_list[i] - self._minimum_training_thermal_parameters[i]
            thermal_parameters_list[i] = thermal_parameters_list[i] / self._maximum_training_thermal_parameters_zeroed[i]
        return thermal_parameters_list

    def predict_A(self, thermal_parameters_list):
        """Predict the value of A (multiplicative correction to photoheating rates) required to give the given thermal
        parameters (T0 in K, gamma and filtering_length in ckpc) at the chosen pivot redshift"""
        A = self._predict_A(*self.map_thermal_parameters_list_to_unit_hypercube(thermal_parameters_list))
        return A

    def predict_B(self, thermal_parameters_list):
        """Predict the value of B (exponent of density-dependent correction to photoheating rates) required to give the
        given thermal parameters (T0 in K, gamma and filtering_length in ckpc) at the chosen pivot redshift"""
        B = self._predict_B(*self.map_thermal_parameters_list_to_unit_hypercube(thermal_parameters_list))
        return B
