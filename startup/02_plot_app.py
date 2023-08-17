import numpy as np
import pandas as pd
# import pkg_resources
from PyQt5 import uic, QtGui, QtCore, QtWidgets
from PyQt5.Qt import QObject

import os
import sys

from lmfit import minimize


path = os.path.dirname(os.path.abspath('__file__'))

import sys


from larch.symboltable import Group

from larch.xafs import xftf

from lmfit import Parameters



import pyqtgraph as pg

from PyQt5.QtCore import pyqtSignal, QThread, QObject

path = os.path.dirname(os.path.abspath('__file__'))



from PyQt5.QtCore import QObject, QThread, pyqtSignal
from time import sleep

class Worker(QObject):
    def __init__(self, k, data, parameters, fit_range, feff_data, kweight=2):
        super(QObject, self).__init__()
        self.k = k
        self.data = data
        self.parameters = parameters
        self.fit_range = fit_range
        self.feff_data = feff_data
        self.kweight = kweight
        self.buffer_parameters = None

    finished = pyqtSignal()
    progress = pyqtSignal(int)

    def calculate_k(self, k, delE0):
        k_new = np.sqrt(k ** 2 - 0.26246592 * delE0)
        return k_new


    def calculate_exafs_lmfit(self, k, parameters, feff_data, fit_range, k_weight=2):

        parameter_dict = parameters.valuesdict()

        if fit_range is None:
            k1 = 0
            k2 = len(k)
        else:
            k1 = np.where(k >= fit_range[0])[0][0]
            k2 = np.where(k >= fit_range[1])[0][0]

        exafs = 0
        for i, key in enumerate(feff_data.keys()):

            feff = self.feff_data[key]['feff']
            R = parameter_dict['r_' + str(i)]
            N = parameter_dict['n_' + str(i)]
            sigma = parameter_dict['ss_' + str(i)]
            delE0 = parameter_dict['e_' + str(i)]
            C3 = parameter_dict['c3_' + str(i)]
            C4 = parameter_dict['c4_' + str(i)]
            S02 = parameter_dict['s02_' + str(i)]

            k_new = self.calculate_k(k, delE0)
            exafs = exafs + (S02 * N *
                             feff['mag[feff]'][k1:k2] / (k_new[k1:k2] * R ** 2) *
                             np.sin(2 * k_new[k1:k2] * R - ((4 / 3) * C3 * k_new[k1:k2] ** 3) + feff['phase[feff]'][
                                                                                                k1:k2] + feff[
                                                                                                             'real[2*phc]'][
                                                                                                         k1:k2]) *
                             np.exp(-2 * k[k1:k2] ** 2 * sigma ** 2) *
                             np.exp((2 / 3) * k[k1:k2] ** 4 * C4) *
                             np.exp(-2 * R / feff['lambda'][k1:k2])) * (k_new[k1:k2] ** k_weight)
        return exafs


    def residual(self, params, x, data, fit_range, feff_data, k_weight=2):

        k1 = np.where(x >= fit_range[0])[0][0]
        k2 = np.where(x >= fit_range[1])[0][0]

        model = self.calculate_exafs_lmfit(x, params, feff_data, fit_range, k_weight)
        return (data[k1:k2] - model)

    def callback(self, params, iter, resid, *args, **kwargs):
        self.buffer_parameters = params

    def run(self):
        """Long-running task."""

        self.out = minimize(self.residual, self.parameters, method='leastsq', iter_cb=self.callback,
                            args=(self.k, self.data, self.fit_range, self.feff_data, 2))

        self.finished.emit()
        print()



class plot_app(*uic.loadUiType(path + '/.ipython/profile_profile_xfit/startup/ui/plot_app.ui')):

    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.chi = None

        self.finished = pyqtSignal()
        self.progress = pyqtSignal(int)
        self.thread = QThread()

        self.create_chi_plot()
        self.create_ft_plot()
        self.set_chi_ft_of_sim()

        self.params = Parameters()

        self.shells = {}

        self.pushButton_load_feff.clicked.connect(self.load_feff_files)
        self.listWidget_feff.setSelectionMode(2)

        self.pushButton_plot.clicked.connect(self.plot_chi_ft_of_sim)

        self.pushButton_make_ft_sim.clicked.connect(self.make_ft_of_sim)




        self.pushButton_load_chi.clicked.connect(self.load_chi)
        self.pushButton_make_ft.clicked.connect(self.make_ft_of_raw)

        self.feff_files = {}
        self.feff_data = {}

        self.cursor = QtGui.QCursor()

        self.fitting_in_progress = False

        # self.pushButton.clicked.connect(self.update_plot)

        windows = ['sine', 'hanning']

        for key in ["", '_sim']:
            getattr(self, 'comboBox_window' + key).addItems(windows)

        self.comboBox_window.addItems(windows)

        # self.spinBox_shells.editingFinished.connect(self.add_n_remove_shells)

        self.shells = {}

        # self.pushButton_load_feff.clicked.connect(self.add_feff_files)
        self.listWidget_feff.setSelectionMode(2)




        self.pushButton_add_shells_to_fit.clicked.connect(self.add_shells_from_feff)
        self.pushButton_fit.clicked.connect(self.perform_exafs_fit)


        self.timer_update_time = QtCore.QTimer(self)
        self.timer_update_time.setInterval(1)
        self.timer_update_time.timeout.connect(self.update_plot)
        self.timer_update_time.start()



    def create_chi_plot(self):

        self.win_chi = pg.GraphicsLayoutWidget()
        self.verticalLayout_chi.addWidget(self.win_chi)
        self.plot_chi = self.win_chi.addPlot()
        self.plot_chi.setTitle('Chi plot')

        self.ref_raw_chi = self.plot_chi.plot()
        self.ref_fit_chi = self.plot_chi.plot()


    def create_ft_plot(self):
        self.win_ft = pg.GraphicsLayoutWidget()
        self.verticalLayout_ft.addWidget(self.win_ft)
        self.plot_ft = self.win_ft.addPlot()
        self.plot_ft.setTitle('FT plot')

        self.ref_raw_ft_mag = self.plot_ft.plot()
        self.ref_raw_ft_img = self.plot_ft.plot()

        self.ref_fit_ft_mag = self.plot_ft.plot()
        self.ref_fit_ft_img = self.plot_ft.plot()



    def set_chi_ft_of_sim(self):
        self.win_chi_sim = pg.GraphicsLayoutWidget()
        self.win_ft_sim = pg.GraphicsLayoutWidget()
        self.horizontalLayout_chi_sim.addWidget(self.win_chi_sim)
        self.horizontalLayout_ft_sim.addWidget(self.win_ft_sim)

        self.plot_chi_sim = self.win_chi_sim.addPlot()
        self.plot_ft_sim = self.win_ft_sim.addPlot()

        self.ref_sim_chi = self.plot_chi_sim.plot()

        self.ref_sim_ft_mag = self.plot_ft_sim.plot()
        self.ref_sim_ft_re = self.plot_ft_sim.plot()

    def load_feff_files(self):
        filespath = QtWidgets.QFileDialog.getOpenFileNames(parent=self,
                                                           directory="/Users/akhiltayal/Library/CloudStorage/GoogleDrive-akhil.tayal.bnl@gmail.com/My Drive/feff_amp_phases/")[0]

        for path in filespath:
            _file_name = os.path.basename(path)
            self.feff_files[_file_name] = path
            _feff, _parameters = self.readfeff(path)
            self.feff_data[_file_name] = {}
            self.feff_data[_file_name]['feff'] = _feff
            self.feff_data[_file_name]['param'] = _parameters

        self.list_widget_items = {}
        self.listWidget_feff.clear()
        for i, key in enumerate(self.feff_files.keys()):
            self.list_widget_items[key] = QtWidgets.QListWidgetItem(key)
            self.listWidget_feff.insertItem(0, self.list_widget_items[key])
            self.list_widget_items[key].setCheckState(0)



    def plot_chi_ft_of_sim(self):
        if getattr(self, 'listWidget_feff').count() > 0:
            status = self.calculate_exafs_from_feff()
            if status:
                self.ref_sim_chi.setData(self.k_sim, np.nan_to_num(self.exafs_sim))
            else:
                print('Nothing is selected')
        else:
            print('Nothing to plot')



    def make_ft_of_sim(self):
        if self.plot_flag:

            test = Group()
            _kweight = self.spinBox_kweight_sim.value()
            win = self.comboBox_window_sim.currentText()
            kmin = self.doubleSpinBox_kmin_sim.value()
            kmax = self.doubleSpinBox_kmax_sim.value()

            xftf(self.k_sim, self.exafs_sim, group=test, kweight=_kweight, window=win, kmin=kmin, kmax=kmax)

            self.ref_sim_chi.setData(self.k_sim, self.exafs_sim * (self.k_sim**_kweight))

            self.ref_sim_ft_mag.setData(test.r, test.chir_mag)
            self.ref_sim_ft_re.setData(test.r, test.chir_re)

        else:
            print('Nothing to Plot')


    def make_ft_of_raw(self):
        test = Group()
        _kweight = self.spinBox_kweight.value()
        win = self.comboBox_window.currentText()
        kmin = self.doubleSpinBox_kmin.value()
        kmax = self.doubleSpinBox_kmax.value()
        xftf(self.chi['k'], self.chi['chi'], group=test, kweight=_kweight, window=win, kmin=kmin, kmax=kmax)


        self.ref_raw_ft_mag.setData(test.r, test.chir_mag)
        self.ref_raw_ft_img.setData(test.r, test.chir_re)

        self.ref_fit_chi.clear()
        self.ref_fit_ft_mag.clear()
        self.ref_fit_ft_img.clear()



    def add_shells_from_feff(self):

        _keys = []

        for i, key in enumerate(self.feff_files.keys()):
            if self.list_widget_items[key].checkState():
                _keys.append(key)

        self.add_n_remove_shells(_keys)
        self.populate_shells_with_default_params(_keys)


    def populate_shells_with_default_params(self, keys):

        for i, key in enumerate(keys):
            if i == 0:
                getattr(self.shells[key], 'doubleSpinBox_' + 'r_value').setValue(self.feff_data[key]['param'][2])
                getattr(self.shells[key], 'doubleSpinBox_' + 'r_max').setValue(self.feff_data[key]['param'][2]*1.05)
                getattr(self.shells[key], 'doubleSpinBox_' + 'r_min').setValue(self.feff_data[key]['param'][2]*0.95)

                getattr(self.shells[key], 'doubleSpinBox_' + 'n_value').setValue(self.feff_data[key]['param'][1])
                getattr(self.shells[key], 'doubleSpinBox_' + 'n_max').setValue(self.feff_data[key]['param'][1]*1.05)
                getattr(self.shells[key], 'doubleSpinBox_' + 'n_min').setValue(self.feff_data[key]['param'][1]*0.95)

                getattr(self.shells[key], 'doubleSpinBox_' + 'ss_value').setValue(0.05)
                getattr(self.shells[key], 'doubleSpinBox_' + 'ss_max').setValue(0.2)
                getattr(self.shells[key], 'doubleSpinBox_' + 'ss_min').setValue(0.01)
                getattr(self.shells[key], 'doubleSpinBox_' + 'c3_value').setValue(0)
                getattr(self.shells[key], 'doubleSpinBox_' + 'c3_max').setValue(0.02)
                getattr(self.shells[key], 'doubleSpinBox_' + 'c3_min').setValue(0)

                getattr(self.shells[key], 'doubleSpinBox_' + 'c4_value').setValue(0)
                getattr(self.shells[key], 'doubleSpinBox_' + 'c4_max').setValue(0.02)
                getattr(self.shells[key], 'doubleSpinBox_' + 'c4_min').setValue(0)

                getattr(self.shells[key], 'doubleSpinBox_' + 'e_value').setValue(0)
                getattr(self.shells[key], 'doubleSpinBox_' + 'e_max').setValue(20)
                getattr(self.shells[key], 'doubleSpinBox_' + 'e_min').setValue(-20)

                getattr(self.shells[key], 'doubleSpinBox_' + 's02_value').setValue(0.8)
                getattr(self.shells[key], 'doubleSpinBox_' + 's02_max').setValue(1.0)
                getattr(self.shells[key], 'doubleSpinBox_' + 's02_min').setValue(0.7)

            else:
                getattr(self.shells[key], 'doubleSpinBox_' + 'r_value').setValue(self.feff_data[key]['param'][2])
                getattr(self.shells[key], 'doubleSpinBox_' + 'r_max').setValue(self.feff_data[key]['param'][2] * 1.05)
                getattr(self.shells[key], 'doubleSpinBox_' + 'r_min').setValue(self.feff_data[key]['param'][2] * 0.95)

                getattr(self.shells[key], 'doubleSpinBox_' + 'n_value').setValue(self.feff_data[key]['param'][1])
                getattr(self.shells[key], 'doubleSpinBox_' + 'n_max').setValue(self.feff_data[key]['param'][1] * 1.05)
                getattr(self.shells[key], 'doubleSpinBox_' + 'n_min').setValue(self.feff_data[key]['param'][1] * 0.95)

                getattr(self.shells[key], 'doubleSpinBox_' + 'ss_value').setValue(0.05)
                getattr(self.shells[key], 'doubleSpinBox_' + 'ss_max').setValue(0.2)
                getattr(self.shells[key], 'doubleSpinBox_' + 'ss_min').setValue(0.01)
                getattr(self.shells[key], 'lineEdit_' + 'ss').setText(f'ss_{i-1:d}')

                getattr(self.shells[key], 'doubleSpinBox_' + 'c3_value').setValue(0)
                getattr(self.shells[key], 'doubleSpinBox_' + 'c3_max').setValue(0.02)
                getattr(self.shells[key], 'doubleSpinBox_' + 'c3_min').setValue(0)

                getattr(self.shells[key], 'doubleSpinBox_' + 'c4_value').setValue(0)
                getattr(self.shells[key], 'doubleSpinBox_' + 'c4_max').setValue(0.02)
                getattr(self.shells[key], 'doubleSpinBox_' + 'c4_min').setValue(0)

                getattr(self.shells[key], 'doubleSpinBox_' + 'e_value').setValue(0)
                getattr(self.shells[key], 'doubleSpinBox_' + 'e_max').setValue(20)
                getattr(self.shells[key], 'doubleSpinBox_' + 'e_min').setValue(-20)
                getattr(self.shells[key], 'lineEdit_' + 'e').setText(f'e_0')

                getattr(self.shells[key], 'doubleSpinBox_' + 's02_value').setValue(0.8)
                getattr(self.shells[key], 'doubleSpinBox_' + 's02_max').setValue(1.0)
                getattr(self.shells[key], 'doubleSpinBox_' + 's02_min').setValue(0.7)
                getattr(self.shells[key], 'lineEdit_' + 's02').setText(f's02_0')


    def load_chi(self):
        self.fitting_in_progress = False
        file = QtWidgets.QFileDialog.getOpenFileName(parent=self)[0]
        col_names = ['k', 'chi', 'chik', 'chik2', 'chik3', 'win', 'energy']
        self.chi = pd.read_csv(file, comment='#', sep='\s+', names=col_names)
        self.plot_raw_chi()
        self.make_ft_of_raw()


    def plot_raw_chi(self):

        _kweight = self.spinBox_kweight.value()
        _chi = self.chi['chi'] *self.chi['k']**_kweight

        self.ref_raw_chi.setData(self.chi['k'], _chi, mkPen='cyan', clear=True)

        self.ref_fit_chi.clear()
        self.ref_fit_ft_mag.clear()
        self.ref_fit_ft_img.clear()




    def perform_exafs_fit(self):
        self.pushButton_fit.setEnabled(False)
        self.params = Parameters()
        self.params.clear()

        self.shells_feff_data = {}


        for i, key in enumerate(self.shells.keys()):

            for par in ['r', 'n', 'ss', 'c3', 'c4', 'e', 's02']:
                _value = getattr(self.shells[key], 'doubleSpinBox_' + par + '_value').value()
                _min = getattr(self.shells[key], 'doubleSpinBox_' + par + '_min').value()
                _max = getattr(self.shells[key], 'doubleSpinBox_' + par + '_max').value()
                _vary = getattr(self.shells[key], 'checkBox_' + par).isChecked()
                _expr = getattr(self.shells[key], 'lineEdit_' + par).text()

                if i == 0 or _expr == "":
                    _expr = None

                self.params.add(f'{par}_{i:d}', value=_value, min=_min, max=_max, vary=_vary, expr=_expr)

        kmin = self.doubleSpinBox_kmin.value()
        kmax = self.doubleSpinBox_kmax.value()

        for key in self.shells.keys():
            self.shells_feff_data[key] = self.feff_data[key]



        self.fit_range = [kmin, kmax]

        self.fitting_in_progress = True

        self.thread = QThread()

        self.worker = Worker(self.chi['k'], self.chi['chik2'], self.params, self.fit_range, self.shells_feff_data)

        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.get_finished_status)
        self.thread.finished.connect(self.thread.deleteLater)
        # self.worker.progress.connect(self.reportProgress)
        self.thread.start()

        print(self.fitting_in_progress)


    def get_finished_status(self):
        print('Fitting finished')
        self.fitting_in_progress = False
        self.pushButton_fit.setEnabled(True)



    def plot_fit_chi_ft(self, parameters=None):

        _parameters = parameters

        exafs_fit = self.calculate_exafs_lmfit(self.chi['k'], _parameters, self.shells.keys(), fit_range=None, k_weight=2)
        _buffer = Group()

        xftf(self.chi['k'], np.nan_to_num(exafs_fit), group=_buffer, kweight=0, window='sine', kmin=3, kmax=14)


        self.ref_fit_chi.setData(self.chi['k'], np.nan_to_num(exafs_fit), pen='yellow')

        self.ref_fit_ft_mag.setData(_buffer.r, _buffer.chir_mag, pen='yellow')
        self.ref_fit_ft_img.setData(_buffer.r, _buffer.chir_re, pen='yellow')


    def update_plot(self):
        if self.fitting_in_progress and (self.worker.buffer_parameters is not None):
            self.plot_fit_chi_ft(parameters=self.worker.buffer_parameters)
            self.populate_shells_with_fit_params(parameters=self.worker.buffer_parameters)



    def populate_shells_with_fit_params(self, parameters=None):

        _parameters = parameters

        _variables = ['r', 'n', 'ss', 'c3', 'c4', 'e', 's02']

        for i, key in enumerate(self.shells.keys()):
            for _var in _variables:
                getattr(self.shells[key], 'doubleSpinBox_' + _var + '_value').setValue(_parameters.valuesdict()[f'{_var}_{i:d}'])

    def calculate_exafs_from_feff(self):
        self.exafs_sim = 0
        self.k_sim = np.arange(0.05, 20, 0.05)
        self.plot_flag = False
        for i, key in enumerate(self.feff_files.keys()):
            if self.list_widget_items[key].checkState():
                self.plot_flag = True
                parameter_dict = {'R': self.feff_data[key]['param'][2],
                                  'N': self.feff_data[key]['param'][1],
                                  'sigma': 0.02,
                                  'DelE0_': 0,
                                  'C3_': 0,
                                  'C4_': 0,
                                  'S02_': 1}

                self.exafs_sim += self.calculate_exafs_wtparam(self.k_sim, parameter_dict, self.feff_data[key]['feff'],
                                                               fit_range=None, k_weight=0)
        return self.plot_flag


    def add_n_remove_shells(self, keys):

        if self.horizontalLayout_param.count() > 0:
            for key in self.shells.copy().keys():
                self.horizontalLayout_param.removeWidget(self.shells[key])
                del self.shells[key]
        else:
            pass

        self.horizontal_spacer = QtWidgets.QSpacerItem(150, 10, QtWidgets.QSizePolicy.Expanding)

        # self.spaceItem = QSpacerItem(150, 10, QSizePolicy.Expanding)

        for i, key in enumerate(keys):
            self.shells[key] = Shell()
            self.horizontalLayout_param.insertWidget(i, self.shells[key])
            self.shells[key].groupBox.setTitle(f'Shell : {i} : ' + key)

    def calculate_k(self, k, delE0):
        k_new = np.sqrt(k ** 2 - 0.26246592 * delE0)
        return k_new

    def readfeff(self, file):
        columns = ['k', 'real[2*phc]', 'mag[feff]', 'phase[feff]', 'red_factor', 'lambda', 'real[p]']

        parameters = []  # parameters nleg, deg, reff, rnrmav(bohr), edge

        with open(file, 'r') as param:
            all_data = param.readlines()
            c = 0
            for i in range(len(all_data)):
                b = all_data[i].find('nleg, deg')
                if b != -1:
                    c = i
                    break
            parameters = []
            for data in all_data[c].split():
                try:
                    parameters.append(float(data))
                except ValueError:
                    parameters.append(data)

            for i in range(c, len(all_data)):
                b = all_data[i].find('mag[feff]')
                if b != -1:
                    c = i
                    break

            feff = pd.read_csv(file, comment='#', skiprows=c + 1, sep='\s+', names=columns)

        return feff, parameters

    def rename_feff_files(self, path, filelist, path_save=None):
        if path_save is None:
            path_save = path
        for file in filelist:
            indx_1 = None
            indx_2 = None
            empty = []
            with open(path + file) as f:
                header = f.readlines()
                name = header[0].split()[0]
                empty.append(name)

                for indx, line in enumerate(header):
                    if 'nleg' in line:
                        data = line.split()

                        _R = data[2]
                        _R_to_include = float(data[2]) * 100
                        _R = str(int(_R_to_include))
                        empty.append(_R)

                        _CN = str(int(float(data[1])))
                        empty.append(_CN)

                    elif 'absorbing' in line:
                        indx_1 = indx
                    elif 'phase' in line:
                        indx_2 = indx


                    elif indx_1 and indx_2:
                        for i in range(indx_1, indx_2):
                            data = header[i].split()[5]
                            empty.append(data)
                        print(empty)
                        break
            rename = ""
            for name in empty:
                rename += name + "_"
            rename = rename[:-1] + ".dat"
            print(f"{rename = }")
            os.rename(os.path.join(path, file), os.path.join(path_save, rename))

    def calculate_exafs_wtparam(self, k, parameter_dict, feff, fit_range, k_weight=2):

        # chi = S02 x N x (f(k)/k x R^2) x sin(2kR + 4/3 x C3 x k^3 + ph(k) + phi_c) x exp(-2k^2sigma^2) * exp(2/3xk^4*C4) * exp(-2R/f(lambda)

        if fit_range is None:
            k1 = 0
            k2 = len(k)
        else:
            k1 = np.where(k >= fit_range[0])[0][0]
            k2 = np.where(k >= fit_range[1])[0][0]

        exafs = 0

        R = parameter_dict['R']
        N = parameter_dict['N']
        sigma = parameter_dict['sigma']
        delE0 = parameter_dict['DelE0_']
        C3 = parameter_dict['C3_']
        C4 = parameter_dict['C4_']
        S02 = parameter_dict['S02_']

        k_new = self.calculate_k(k, delE0)
        exafs = exafs + (S02 * N *
                         feff['mag[feff]'][k1:k2] / (k_new[k1:k2] * R ** 2) *
                         np.sin(2 * k_new[k1:k2] * R + ((4 / 3) * C3 * k_new[k1:k2] ** 3) + feff['phase[feff]'][k1:k2] +
                                feff['real[2*phc]'][k1:k2]) *
                         np.exp(-2 * k[k1:k2] ** 2 * sigma ** 2) *
                         np.exp((2 / 3) * k[k1:k2] ** 4 * C4) *
                         np.exp(-2 * R / feff['lambda'][k1:k2])) * (k_new[k1:k2] ** k_weight)
        return exafs

    def calculate_exafs_lmfit(self, k, parameters, feff_keys, fit_range, k_weight=2):

        parameter_dict = parameters.valuesdict()

        if fit_range is None:
            k1 = 0
            k2 = len(k)
        else:
            k1 = np.where(k >= fit_range[0])[0][0]
            k2 = np.where(k >= fit_range[1])[0][0]

        exafs = 0
        for i, key in enumerate(feff_keys):

            feff = self.feff_data[key]['feff']
            R = parameter_dict['r_' + str(i)]
            N = parameter_dict['n_' + str(i)]
            sigma = parameter_dict['ss_' + str(i)]
            delE0 = parameter_dict['e_' + str(i)]
            C3 = parameter_dict['c3_' + str(i)]
            C4 = parameter_dict['c4_' + str(i)]
            S02 = parameter_dict['s02_' + str(i)]

            k_new = self.calculate_k(k, delE0)
            exafs = exafs + (S02 * N *
                             feff['mag[feff]'][k1:k2] / (k_new[k1:k2] * R ** 2) *
                             np.sin(2 * k_new[k1:k2] * R - ((4 / 3) * C3 * k_new[k1:k2] ** 3) + feff['phase[feff]'][
                                                                                                k1:k2] + feff[
                                                                                                             'real[2*phc]'][
                                                                                                         k1:k2]) *
                             np.exp(-2 * k[k1:k2] ** 2 * sigma ** 2) *
                             np.exp((2 / 3) * k[k1:k2] ** 4 * C4) *
                             np.exp(-2 * R / feff['lambda'][k1:k2])) * (k_new[k1:k2] ** k_weight)
        return exafs

    def residual(self, params, x, data, fit_range, feff_keys, k_weight=2):

        k1 = np.where(x >= fit_range[0])[0][0]
        k2 = np.where(x >= fit_range[1])[0][0]

        model = self.calculate_exafs_lmfit(x, params, feff_keys, fit_range, k_weight)
        return (data[k1:k2] - model)





    def save_data_files(self, file_name, out, path_to_save, feff_files, data_file, k, chi, k_weight=2,
                        k_min=3, k_max=14):
        with open(path_to_save + file_name + ".INP", 'w+') as fnm:
            fnm.writelines(["Data File = ", data_file, "\n\n",
                            "Method = ", out.method, "\n",
                            "chi-squre = ", str(out.chisqr), "\n",
                            "Reduced chi-square = ", str(out.redchi), "\n",
                            "Akaike info crit. = ", str(out.aic), "\n",
                            "Bayesian info crit. = ", str(out.bic), "\n",
                            "No. of variables = ", str(out.nvarys), "\n\n",
                            ])
            parameters = out.params.valuesdict()

            for i in range(len(feff_files)):
                fnm.writelines(["Shell ", str(i), "\n"])
                if out.params["R" + str(i)].stderr is not None:

                    fnm.writelines(["feff file = ", feff_files[i], "\n\n",
                                    f"R = {parameters['R' + str(i)]:2.3f}+/- {out.params['R' + str(i)].stderr:2.3f}\n",
                                    f"N = {parameters['N' + str(i)]:2.3f}+/- {out.params['N' + str(i)].stderr:2.2f}\n",
                                    f"sigma = {parameters['sigma' + str(i)]:2.3f}+/- {out.params['sigma' + str(i)].stderr:2.3f}\n",
                                    f"C3 = {parameters['C3_' + str(i)]:2.5f}+/- {out.params['C3_' + str(i)].stderr:2.5f}\n"
                                    f"C4 = {parameters['C4_' + str(i)]:2.6f}+/- {out.params['C4_' + str(i)].stderr:2.6f}\n"
                                    f"DelE0 = {parameters['DelE0_' + str(i)]:2.3f}+/- {out.params['DelE0_' + str(i)].stderr:2.2f}\n",
                                    f"S02 = {parameters['S02_' + str(i)]:2.3f}+/- {out.params['S02_' + str(i)].stderr:2.3f}\n\n"])

                else:
                    fnm.writelines(["feff file = ", feff_files[i], "\n\n",
                                    f"R = {parameters['R' + str(i)]:2.3f}+/- {0}\n",
                                    f"N = {parameters['N' + str(i)]:2.3f}+/- {0}\n",
                                    f"sigma = {parameters['sigma' + str(i)]:2.3f}+/- {0}\n",
                                    f"C3 = {parameters['C3_' + str(i)]:2.5f}+/- {0}\n"
                                    f"C4 = {parameters['C4_' + str(i)]:2.6f}+/- {0}\n"
                                    f"DelE0 = {parameters['DelE0_' + str(i)]:2.3f}+/- {0}\n",
                                    f"S02 = {parameters['S02_' + str(i)]:2.3f}+/- {0}\n\n"])

        parameters = out.params.valuesdict()

        arr = [np.array(k)]
        arr.append(np.nan_to_num(chi))

        arr_chir = []
        fit = Group()
        xftf(k, np.nan_to_num(chi * k ** k_weight), group=fit, window='sine', kmin=k_min, kmax=k_max)

        arr_chir.append(np.array(fit.r))

        arr_chir.append(np.array(fit.chir_mag))
        arr_chir.append(np.array(fit.chir_re))

        buff_chir = fit.chir_mag

        exafs = self.calculate_exafs_lmfit(k, out.params, feff_files, None, 0)

        del (fit)
        fit = Group()
        xftf(k, np.nan_to_num(exafs * k ** k_weight), group=fit, window='sine', kmin=k_min, kmax=k_max)
        arr_chir.append(np.array(fit.chir_mag))
        arr_chir.append(np.array(fit.chir_re))

        arr_chir.append(np.subtract(buff_chir, fit.chir_mag))

        arr.append(np.nan_to_num(exafs))

        arr.append(np.nan_to_num(np.subtract(chi, exafs)))

        indx = []
        # np.savetxt(path2+"\\"+file_name+".FIT", np.row_stack(dat.k))
        for i in range(len(feff_files)):
            par_dict = {key[:-1]: out.params.valuesdict()[key] for key in out.params.valuesdict() if key[-1] == str(i)}
            exafs = self.calculate_exafs_wtparam(k, par_dict, feff_files[i], None, k_weight=0)

            arr.append(np.nan_to_num(exafs))

            del (fit)
            fit = Group()
            xftf(k, np.nan_to_num(exafs * k ** k_weight), group=fit, window='sine', kmin=k_min, kmax=k_max)

            arr_chir.append(np.array(fit.chir_mag))

            indx.append(i)
        tup = tuple(arr)

        tup_r = tuple(arr_chir)

        wave = ''
        for i in indx:
            wave += 'wave' + str(i + 1) + " "

        header_chi = "k chi_data chi_fit chi_diff " + wave
        header_chir = "R chir_data chir_datar chir_fit chir_fitr chir_diff " + wave

        np.savetxt(path_to_save + file_name + ".FIT", np.column_stack(tup), header=header_chi)
        np.savetxt(path_to_save + file_name + ".FTM", np.column_stack(tup_r), header=header_chir)

if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    window = plot_app()
    window.show()
    sys.exit(app.exec_())