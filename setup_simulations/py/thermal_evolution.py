import os
import numpy as np


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
