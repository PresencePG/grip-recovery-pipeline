import pandas as pd
import numpy as np
import re
from datetime import datetime as dt
from dateutil import parser
import cmath
import json
import os


def str_to_isoformat_str(s):
    dateutil_time = parser.parse(s)
    isoformat_str = dateutil_time.isoformat()
    return isoformat_str


def str_to_isoformat_datetime(s):
    return parser.parse(str_to_isoformat_str(s))


def on_init(t):
    # # Set initial conditions and initiate global parameters
    # global ES_2D, targets, ESCs, starttime, stoptime, nc
    #
    # ES_2D = True
    #
    # # Read ESCs,targets and time data from json file
    # with open('config.json', 'r') as json_file:
    #     data = json.load(json_file)
    #
    # starttime_dt = data['starttime']
    # starttime = int((dt.strptime(starttime_dt, '%Y-%m-%dT%H:%M:%S-00:00') - dt(1970,1,1)).total_seconds())
    #
    # stoptime_dt = data['stoptime']
    # stoptime = int((dt.strptime(stoptime_dt, '%Y-%m-%dT%H:%M:%S-00:00') - dt(1970,1,1)).total_seconds())
    #
    # ESCs = data['controllers']
    # targets = data['optimization_targets']
    #
    # nc = len(ESCs) # Number of connected ES controllers

    return True


def on_commit(t):
    # This computer time will be compared with the computer time after the iteration to measure the computation time of an iteration.
    time_before_step = dt.now()

    # Check if iteration_count exists already. If not, initialize parameter and simulation
    if 'iteration_count' not in globals():

        # Initialize iteration_count
        global iteration_count
        iteration_count = 0

        # Set initial conditions and initiate global parameters
        global ES_2D, targets, ESCs, starttime, stoptime, nc, totaltime, t0

        ES_2D = True

        # Read ESCs,targets and time data from json file
        with open('config.json', 'r') as json_file:
            data = json.load(json_file)

        starttime_dt = data['starttime']
        starttime = int((str_to_isoformat_datetime(starttime_dt) - str_to_isoformat_datetime(
            '1970-01-01T00:00:00+00:00')).total_seconds())

        stoptime_dt = data['stoptime']
        stoptime = int((str_to_isoformat_datetime(stoptime_dt) - str_to_isoformat_datetime(
            '1970-01-01T00:00:00+00:00')).total_seconds())

        totaltime = int(
            (str_to_isoformat_datetime(stoptime_dt) - str_to_isoformat_datetime(starttime_dt)).total_seconds())

        #print('STARTTIME')
        #print((str_to_isoformat_datetime(starttime_dt)))
        #print(starttime)

        #print('STOPTIME')
        #print((str_to_isoformat_datetime(stoptime_dt)))
        #print(stoptime)

        t0 = t

        #print('TO')
        #print(t0)

        ESCs = data['controllers']
        targets = data['optimization_targets']

        nc = len(ESCs)  # Number of connected ES controllers

        global results, results_all, columns_results
        results_all = []

        ''' Create all columns for the results dataframe as a list. Columns here need to exactly match the order of the results list created later in the code '''
        columns_results = ['UNIX']
        columns_results.append('# timestamp')
        for i in range(len(targets)):
            columns_results.append('{}:{}'.format(targets[i]['node_1_name'], 'measured_voltage_A[V]'))
            columns_results.append('{}:{}'.format(targets[i]['node_1_name'], 'measured_voltage_B[V]'))
            columns_results.append('{}:{}'.format(targets[i]['node_1_name'], 'measured_voltage_C[V]'))
            columns_results.append('{}:{}'.format(targets[i]['node_2_name'], 'measured_voltage_A[V]'))
            columns_results.append('{}:{}'.format(targets[i]['node_2_name'], 'measured_voltage_B[V]'))
            columns_results.append('{}:{}'.format(targets[i]['node_2_name'], 'measured_voltage_C[V]'))
            # for phase in targets[i]['phases']:
            for phase in ['A', 'B', 'C']:
                columns_results.append('angle_{}_{}'.format(targets[i]['node_1_name'], phase))
                columns_results.append('angle_{}_{}'.format(targets[i]['node_2_name'], phase))
                columns_results.append(
                    'diff_angle_{}_{}_{}'.format(targets[i]['node_1_name'], targets[i]['node_2_name'], phase))
                columns_results.append(
                    'V_diff_{}_{}_{}'.format(targets[i]['node_1_name'], targets[i]['node_2_name'], phase))
            # columns_results.append('line_power_{}_{}'.format('702','705')) # Measured line power is only for debugging. Can't be used in the simulations, because target nodes might not be neighboring nodes, so they could be connected by multiple lines in a row.
            for ESC in ESCs:
                columns_results.append(
                    '{}:measured_power_{}.real[kW]'.format(ESC['node_name'], ESC['controller_phase']))
                if ES_2D:
                    columns_results.append(
                        '{}:measured_power_{}.imag[kVAr]'.format(ESC['node_name'], ESC['controller_phase']))
            for i in range(nc):
                columns_results.append('P_ES_{}'.format(i + 1))
                if ES_2D:
                    columns_results.append('Q_ES_{}'.format(i + 1))
        columns_results.append('objective_function_value')

        global T, time, computertime_0
        computertime_0 = dt.now()
        T = 1.0
        time = np.arange(0, stoptime - starttime + T, T)

        # objective function
        global J
        J = np.zeros(len(time))

        ''' SETUP REGULAR ES CONTROLLERS '''
        global rho, eps, sigma, xi, xiavg, uhat, u
        # setup vectors for ES controller internal states
        rho = np.zeros((nc, len(time)))  # value after highpass filter
        eps = np.zeros((nc, len(time)))  # LPF of objective
        sigma = np.zeros((nc, len(time)))  # value after demodulation
        xi = np.zeros((nc, len(time)))  # value after lowpass filter
        xiavg = np.zeros((nc, len(time)))  # averaged value after lowpass filter
        uhat = np.zeros((nc, len(time)))  # setpoint - value after integration
        u = np.zeros((nc, len(time)))  # controller output - control outputs of ES control (position)

        # controller parameters
        global f, wes, a, hpf, whpf, lpf, wlpf, kint
        f = np.zeros(nc)
        a = np.zeros(nc)
        kint = np.zeros(nc)
        for i in range(len(ESCs)):
            f[i] = ESCs[i]['controller_frequency']  # frequency of probe [Hz]
            a[i] = ESCs[i]['real_power_probing_amplitude'] * 1000  # probe amplitude
            kint[i] = ESCs[i]['controller_speed'] * (-1)  # gain on the integrator
        wes = 2 * np.pi * f  # angular frequency of probe
        hpf = np.array(f / 10)  # highpass filter critical frequency [Hz]
        whpf = 2 * np.pi * hpf  # angular frequency of highpass filter
        lpf = np.array(f / 10)  # lowpass filter critical frequency [Hz]
        wlpf = 2 * np.pi * lpf  # angular frequency of lowwpass filter

        # controller initial conditions
        uhat[:, 0] = np.zeros(nc)
        u[:, 0] = uhat[:, 0]
        eps[:, 0] = 0
        rho[:, 0] = 0
        sigma[:, 0] = 0
        xi[:, 0] = 0

        ''' SETUP 2D ES CONTROLLERS '''
        if ES_2D:
            global rho2D, eps2D, sigma2D, xi2D, xiavg2D, uhat2D, u2D
            # setup vectors for ES controller internal states
            rho2D = np.zeros((nc, len(time)))  # value after highpass filter
            eps2D = np.zeros((nc, len(time)))  # LPF of objective
            sigma2D = np.zeros((nc, len(time)))  # value after demodulation
            xi2D = np.zeros((nc, len(time)))  # value after lowpass filter
            xiavg2D = np.zeros((nc, len(time)))  # averaged value after lowpass filter
            uhat2D = np.zeros((nc, len(time)))  # setpoint - value after integration
            u2D = np.zeros((nc, len(time)))  # controller output - control outputs of ES control (position)

            # controller 2D parameters
            global wes2D, a2D, hpf2D, whpf2D, lpf2D, wlpf2D, kint2D
            a2D = np.zeros(nc)
            kint2D = np.zeros(nc)
            for i in range(len(ESCs)):
                a2D[i] = ESCs[i]['reactive_power_probing_amplitude'] * 1000  # probe amplitude
                kint2D[i] = ESCs[i]['controller_speed'] * (-1)  # gain on the integrator
            wes2D = 2 * np.pi * f  # angular frequency of probe
            hpf2D = np.array(f / 10)  # highpass filter critical frequency [Hz]
            whpf2D = 2 * np.pi * hpf  # angular frequency of highpass filter
            lpf2D = np.array(f / 10)  # lowpass filter critical frequency [Hz]
            wlpf2D = 2 * np.pi * lpf  # angular frequency of lowwpass filter

            # controller initial conditions
            uhat2D[:, 0] = np.zeros(nc)
            u2D[:, 0] = uhat2D[:, 0]
            eps2D[:, 0] = 0
            rho2D[:, 0] = 0
            sigma2D[:, 0] = 0
            xi2D[:, 0] = 0

    # Set kt to current value
    # kt = t - starttime
    kt = t - t0

    # Print information for debugging about timesteps
    # print(dt.fromtimestamp(int(t)).strftime('%Y-%m-%dT%H:%M:%S-00:00'))
    # print(starttime)
    # print(t)

    # compute objective function values
    # if t >= starttime:
    if t - t0 >= 0:
        grid_outputs = {}

        '''Get values from gridlabd for timestep and save in dict'''
        for target in targets:
            grid_outputs['V_{}_{}'.format(target['node_1_name'], 'A')] = complex(
                re.sub("V", "", gridlabd.get_value('{}'.format(target['node_1_name']), 'voltage_A')))
            grid_outputs['V_{}_{}'.format(target['node_1_name'], 'B')] = complex(
                re.sub("V", "", gridlabd.get_value('{}'.format(target['node_1_name']), 'voltage_B')))
            grid_outputs['V_{}_{}'.format(target['node_1_name'], 'C')] = complex(
                re.sub("V", "", gridlabd.get_value('{}'.format(target['node_1_name']), 'voltage_C')))
            grid_outputs['V_{}_{}'.format(target['node_2_name'], 'A')] = complex(
                re.sub("V", "", gridlabd.get_value('{}'.format(target['node_2_name']), 'voltage_A')))
            grid_outputs['V_{}_{}'.format(target['node_2_name'], 'B')] = complex(
                re.sub("V", "", gridlabd.get_value('{}'.format(target['node_2_name']), 'voltage_B')))
            grid_outputs['V_{}_{}'.format(target['node_2_name'], 'C')] = complex(
                re.sub("V", "", gridlabd.get_value('{}'.format(target['node_2_name']), 'voltage_C')))

            grid_outputs['S_{}_{}'.format(target['node_1_name'], 'A')] = abs(
                grid_outputs['V_{}_{}'.format(target['node_1_name'], 'A')])
            grid_outputs['S_{}_{}'.format(target['node_1_name'], 'B')] = abs(
                grid_outputs['V_{}_{}'.format(target['node_1_name'], 'B')])
            grid_outputs['S_{}_{}'.format(target['node_1_name'], 'C')] = abs(
                grid_outputs['V_{}_{}'.format(target['node_1_name'], 'C')])
            grid_outputs['S_{}_{}'.format(target['node_2_name'], 'A')] = abs(
                grid_outputs['V_{}_{}'.format(target['node_2_name'], 'A')])
            grid_outputs['S_{}_{}'.format(target['node_2_name'], 'B')] = abs(
                grid_outputs['V_{}_{}'.format(target['node_2_name'], 'B')])
            grid_outputs['S_{}_{}'.format(target['node_2_name'], 'C')] = abs(
                grid_outputs['V_{}_{}'.format(target['node_2_name'], 'C')])

            grid_outputs['angle_{}_{}'.format(target['node_1_name'], 'A')] = cmath.phase(
                grid_outputs['V_{}_{}'.format(target['node_1_name'], 'A')])
            grid_outputs['angle_{}_{}'.format(target['node_1_name'], 'B')] = cmath.phase(
                grid_outputs['V_{}_{}'.format(target['node_1_name'], 'B')])
            grid_outputs['angle_{}_{}'.format(target['node_1_name'], 'C')] = cmath.phase(
                grid_outputs['V_{}_{}'.format(target['node_1_name'], 'C')])
            grid_outputs['angle_{}_{}'.format(target['node_2_name'], 'A')] = cmath.phase(
                grid_outputs['V_{}_{}'.format(target['node_2_name'], 'A')])
            grid_outputs['angle_{}_{}'.format(target['node_2_name'], 'B')] = cmath.phase(
                grid_outputs['V_{}_{}'.format(target['node_2_name'], 'B')])
            grid_outputs['angle_{}_{}'.format(target['node_2_name'], 'C')] = cmath.phase(
                grid_outputs['V_{}_{}'.format(target['node_2_name'], 'C')])

            grid_outputs['V_diff_{}_{}_{}'.format(target['node_1_name'], target['node_2_name'], 'A')] = grid_outputs[
                                                                                                            'V_{}_{}'.format(
                                                                                                                target[
                                                                                                                    'node_1_name'],
                                                                                                                'A')] - \
                                                                                                        grid_outputs[
                                                                                                            'V_{}_{}'.format(
                                                                                                                target[
                                                                                                                    'node_2_name'],
                                                                                                                'A')]
            grid_outputs['V_diff_{}_{}_{}'.format(target['node_1_name'], target['node_2_name'], 'B')] = grid_outputs[
                                                                                                            'V_{}_{}'.format(
                                                                                                                target[
                                                                                                                    'node_1_name'],
                                                                                                                'B')] - \
                                                                                                        grid_outputs[
                                                                                                            'V_{}_{}'.format(
                                                                                                                target[
                                                                                                                    'node_2_name'],
                                                                                                                'B')]
            grid_outputs['V_diff_{}_{}_{}'.format(target['node_1_name'], target['node_2_name'], 'C')] = grid_outputs[
                                                                                                            'V_{}_{}'.format(
                                                                                                                target[
                                                                                                                    'node_1_name'],
                                                                                                                'C')] - \
                                                                                                        grid_outputs[
                                                                                                            'V_{}_{}'.format(
                                                                                                                target[
                                                                                                                    'node_2_name'],
                                                                                                                'C')]

            grid_outputs['angle_diff_{}_{}_{}'.format(target['node_1_name'], target['node_2_name'], 'A')] = \
            grid_outputs['angle_{}_{}'.format(target['node_1_name'], 'A')] - grid_outputs[
                'angle_{}_{}'.format(target['node_2_name'], 'A')]
            grid_outputs['angle_diff_{}_{}_{}'.format(target['node_1_name'], target['node_2_name'], 'B')] = \
            grid_outputs['angle_{}_{}'.format(target['node_1_name'], 'B')] - grid_outputs[
                'angle_{}_{}'.format(target['node_2_name'], 'B')]
            grid_outputs['angle_diff_{}_{}_{}'.format(target['node_1_name'], target['node_2_name'], 'C')] = \
            grid_outputs['angle_{}_{}'.format(target['node_1_name'], 'C')] - grid_outputs[
                'angle_{}_{}'.format(target['node_2_name'], 'C')]

            # Fixed line power, not flexible to chose by the user
            # grid_outputs['line_power_702_705'] = gridlabd.get_value('line_702_705', 'power_in')

        ''' *** OBJECTIVE FUNCTION *** '''
        Jk_phase = 0
        Jk_magnitude = 0

        for target in targets:
            # for phase in target['phases']:
            for phase in ['A', 'B', 'C']:
                Jk_phase += (grid_outputs[
                                 'angle_diff_{}_{}_{}'.format(target['node_1_name'], target['node_2_name'], phase)] *
                             target[
                                 'weighting_factor'] * 10000) ** 2  # The multiplication with 10000 is because the phase differences are on a much smaller scale.

        for target in targets:
            # for phase in target['phases']:
            for phase in ['A', 'B', 'C']:
                Jk_magnitude += (grid_outputs['V_diff_{}_{}_{}'.format(target['node_1_name'], target['node_2_name'],
                                                                       phase)].real * target['weighting_factor']) ** 2

        Jk = Jk_phase + Jk_magnitude
        J[kt] = Jk

        #print('Timestep: ', kt, '   Objective Function: ', Jk)
        #print('t - t0: ', t - t0, '  totaltime: ', totaltime)

    if kt >= 1:
        # highpass filter
        rho[:, kt] = (1 - whpf * T) * rho[:, kt - 1] + J[kt] - J[kt - 1]

        # objective function error
        eps[:, kt] = J[kt] - rho[:, kt]

        # demodulate
        sigma[:, kt] = 2 / a * np.cos(wes * time[kt - 1]) * rho[:, kt]

        # lowpass filter demodulated values
        xi[:, kt] = (1 - wlpf * T) * xi[:, kt - 1] + wlpf * T * sigma[:, kt - 1]

        # integrate to obtain setpoint
        uhat[:, kt] = uhat[:, kt - 1] + kint * T * xi[:, kt - 1]

        # Read solar generation from gridlabD for every ESC and limit setpoints to generation constraints
        solar_gen = {}
        for i in range(len(ESCs)):

            # Read solar generation from gridlabD
            solar_gen['{}'.format(i)] = complex(
                re.sub("VA", "", gridlabd.get_value('ES_{}_solar'.format(i + 1), 'VA_Out'))).real

            # Limit setpoint to physical limits
            clamp = lambda n, minn, maxn: max(min(maxn, n), minn)
            uhat[i, kt] = clamp(uhat[i, kt], ESCs[i]['real_power_probing_amplitude'] * 1000,
                                solar_gen[str(i)] - ESCs[i]['real_power_probing_amplitude'] * 1000)

            # Check if enough solar generation exist for full probing amplitude. If True, add probe to setpoint, if False, deactivate ES
            if solar_gen[str(i)] >= ESCs[i]['real_power_probing_amplitude'] * 1000:
                # add probe to setpoint
                u[i, kt] = uhat[i, kt] + a[i] * np.cos(wes[i] * time[kt])
            else:
                # If maximum probe amplitude cannot be utilized, deactivate ES outputs entirely
                u[i, kt] = 0
                uhat[i, kt] = 0

        # print('old solar_gen to be replaced with dictionary')
        # solar_gen = complex(re.sub("VA","",gridlabd.get_value('ES_1_solar','VA_Out'))).real
        #
        # print('needs to be changed to individual controller settings')
        # # limit setpoint
        # uhat = np.clip(uhat, a_min = ESCs[0]['real_power_probing_amplitude'], a_max = 0 + solar_gen[i] - ESCs[0]['real_power_probing_amplitude'])
        #

        for i in range(nc):
            gridlabd.set_value('ES_{}_inv'.format(i + 1), 'P_Out', str("%.6f" % u[i, kt]) + ' W')

        if ES_2D:
            # highpass filter
            rho2D[:, kt] = (1 - whpf2D * T) * rho2D[:, kt - 1] + J[kt] - J[kt - 1]

            # objective function error
            eps2D[:, kt] = J[kt] - rho2D[:, kt]

            # demodulate
            sigma2D[:, kt] = 2 / a2D * np.sin(wes2D * time[kt - 1]) * rho2D[:, kt]

            # lowpass filter demodulated values
            xi2D[:, kt] = (1 - wlpf2D * T) * xi2D[:, kt - 1] + wlpf2D * T * sigma2D[:, kt - 1]

            # integrate to obtain setpoint
            uhat2D[:, kt] = uhat2D[:, kt - 1] + kint2D * T * xi2D[:, kt - 1]

            # limit setpoint to reactive power capacity limit
            uhat2D = np.clip(uhat2D, a_min=-(1000 * ESCs[0]['apparent_power_capacity_limit']) + ESCs[0][
                'reactive_power_probing_amplitude'] * 1000,
                             a_max=1000 * ESCs[0]['reactive_power_capacity_limit'] - ESCs[0][
                                 'reactive_power_probing_amplitude'] * 1000)

            # add probe to setpoint
            u2D[:, kt] = uhat2D[:, kt] + a2D * np.sin(wes2D * time[kt])

            for i in range(nc):
                gridlabd.set_value('ES_{}_inv'.format(i + 1), 'Q_Out', str("%.6f" % u2D[i, kt]) + ' VAr')

        now = dt.now()

        # Create results list for current timestep
        results_step = [t]  # Add Unix time to results
        results_step.append(dt.fromtimestamp(int(t - 1)).strftime(
            '%Y-%m-%dT%H:%M:%S-00:00'))  # Add time in datetime format matching gridlabD timeformat.
        for target in targets:
            # Add to results: Target voltages as complex values
            results_step.append(re.sub("V", "", gridlabd.get_value('{}'.format(target['node_1_name']), 'voltage_A')))
            results_step.append(re.sub("V", "", gridlabd.get_value('{}'.format(target['node_1_name']), 'voltage_B')))
            results_step.append(re.sub("V", "", gridlabd.get_value('{}'.format(target['node_1_name']), 'voltage_C')))
            results_step.append(re.sub("V", "", gridlabd.get_value('{}'.format(target['node_2_name']), 'voltage_A')))
            results_step.append(re.sub("V", "", gridlabd.get_value('{}'.format(target['node_2_name']), 'voltage_B')))
            results_step.append(re.sub("V", "", gridlabd.get_value('{}'.format(target['node_2_name']), 'voltage_C')))
            # for phase in target['phases']:
            for phase in ['A', 'B', 'C']:
                # Add to results: Target node voltage phase angles and voltage phase angle and voltage magnitude differences between target nodes.
                results_step.append(grid_outputs['angle_{}_{}'.format(target['node_1_name'], phase)])
                results_step.append(grid_outputs['angle_{}_{}'.format(target['node_2_name'], phase)])
                results_step.append(
                    grid_outputs['angle_diff_{}_{}_{}'.format(target['node_1_name'], target['node_2_name'], phase)])
                results_step.append(
                    grid_outputs['V_diff_{}_{}_{}'.format(target['node_1_name'], target['node_2_name'], phase)].real)

            # results_step.append(grid_outputs['line_power_702_705']) # Add line power for comparison reason for debugging - not used in GRIP platform.

            # Read inverter outputs
            for i in range(len(ESCs)):
                i = str(i + 1)
                results_step.append(
                    float(re.sub("W", "", gridlabd.get_value('ES_{}_meter'.format(i), 'measured_real_power'))) * (
                        -1) / 1000)
                if ES_2D:
                    results_step.append(float(
                        re.sub("VAr", "", gridlabd.get_value('ES_{}_meter'.format(i), 'measured_reactive_power'))) * (
                                            -1) / 1000)

            for i in range(nc):
                results_step.append(u[i, kt])
                if ES_2D:
                    results_step.append(u2D[i, kt])

        results_step.append(Jk)

        # Append current timestep results to list with results of whole simulation
        results_all.append(results_step)

        # print('Append results:',dt.now() - now)
        # print('Time last timestep:',dt.now() - time_before_step)
        # print('SOLAR POWER GENERATED:', complex(re.sub("VA","",gridlabd.get_value('ES_1_solar','VA_Out'))).real)
        # print('SOLAR OBJECT:', gridlabd.get_object('ES_1_solar'))
        # print('Voltage Microgrid:', gridlabd.get_object('node_705'))
        # print('Voltage 705 A: ', gridlabd.get_value('node_705','voltage_A'))

    # if t == stoptime:
    if t - t0 == totaltime:
        print('columns_results', columns_results)
        print('results_all', results_all)
        # print(results_step)

        results = pd.DataFrame(results_all, columns=columns_results)
        results = results.set_index(results['# timestamp'])
        results = results.drop(['# timestamp'], axis=1)

        print(results)
        # print('Time for whole simulation: ', dt.now() - computertime_0)
        with open('results.csv', mode='w', newline='') as csvfile:
            results.to_csv(csvfile)

    iteration_count += 1  # Increase iteration counter by 1 after finishing iteration
    return True


####### Good to know things below #######
'''
gridlabd.get_object(obj_name) --> this gives you all the information about any kind of objects
"%.6f" % power_b_825 --> You have to have a string of a number with 6 digits after the decimal for gridlabd to work.
'''
