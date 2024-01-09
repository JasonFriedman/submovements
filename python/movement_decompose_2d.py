""" module movement_decompose_2d
A module to decompose a 2d movement to multiple submovements
Contain the following functions:
    load_data            - read a folder full of csv files, 
                           collect movements position, velocities & recorded time
    plot_position        - plot movement position in time
    plot_velocity        - plot movement velocity in time
    decompose_2D         - estimate submovements parameters from movement
    plot_submovements_2D - plot the expected velocities from submovement group

The python code is based on the Matlab code, and was mostly written during a Hackathon:
https://github.com/Liordemarcas/submovements

by
Omer Ophir:             https://github.com/omerophir
Omri FK:                https://github.com/OmriFK
Shay Eylon:             https://github.com/ShayEylon
Lior de Marcas (LdM):   https://github.com/Liordemarcas

Updated by Jason Friedman https://github.com/JasonFriedman
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import filtfilt, butter
from scipy.optimize import minimize

def load_data(dir_name):
    """
    Loads data from CSV files in the specified directory.

    Parameters:
        dir_name (str): Directory path containing the CSV files.
                        In each file expect:
                            Cols 0 & 1 to be X-Y position
                            Col 3 to be pen-pressure (only colecting data when it is >0)
                            Col 4 to be time
    Returns:
        position_filtered (list): List of position data arrays after filtering.
        velocity (list): List of velocity data arrays.
        time (list): List of time data arrays.

    Raises:
        ValueError: If the directory does not contain any CSV files.

    """

    # Get a list of files in the directory
    files = os.listdir(dir_name)
    # Filter only the CSV files
    csv_files = [f for f in files if f.endswith('.csv')]
    # Raise an error if no CSV files are found
    if not csv_files:
        raise ValueError('Must specify a directory to load the csv files from')
    # Extract block and trial information from file names
    blocks = []
    trials = []
    file_names = []
    for file_name in csv_files:
        file_names.append(file_name)
        match = re.search(r'tb_.*block(\d*)_trial(\d*).csv', file_name) #checking for correct file name
        block = int(match.group(1))
        trial = int(match.group(2))
        blocks.append(block)
        trials.append(trial)

    # We have lists of blocks and trials and looks for max to see how much blocks and trials we have in this folder
    max_block = max(blocks)
    max_trial = max(trials)
    position_filtered = []
    velocity = []
    time = []

    # Process data for each block and trial
    for block in range(1, max_block + 1):
        for trial in range(1, max_trial + 1):
            trial_index = [i for i, (_block, _trial) in enumerate(zip(blocks, trials))
                           if _block == block and _trial == trial]
            if not trial_index:
                continue

            data = np.loadtxt(os.path.join(dir_name, csv_files[trial_index[0]]), delimiter=',')
            pressure = data[:, 3]
            position = data[pressure > 0, :2] / 1000
            _time = data[pressure > 0, 4] / 1000  # seconds
            _time = _time - _time[0]
            dt = np.median(np.diff(_time))
            b, a = butter(2, 5 / ((1 / dt) / 2))
            _position_filtered = filtfilt(b, a, position, axis=0)
            _velocity = np.vstack([[0, 0], np.diff(_position_filtered, axis=0) / dt])

            #organizing the data in correct variables for future functions
            time.append(_time)
            position_filtered.append(_position_filtered)
            velocity.append(_velocity)

    return position_filtered, velocity, time

def plot_position(position, time, plot_type=1):
    """
    Plot the position data against time.

    Parameters:
        position (list): List of position data arrays.
        time (list): List of time data arrays.
        plot_type (int, optional): Type of plot to generate. Default is 1.
            - plot_type = 1: x vs y position
            - plot_type = 2: Time vs. x & y position.

    Raises:
        ValueError: If the plot_type is unknown.

    """
    if plot_type not in [1, 2]:
        raise ValueError('Unknown plot type')

    num_positions = len(position)
    cols = int(np.ceil(np.sqrt(num_positions)))
    rows = int(np.ceil(num_positions / cols))
    #create demo figure with num subplots correlated with num of trials in data
    _ , axs = plt.subplots(rows, cols, figsize=(15, 15))
    axs_flat = axs.flatten()
    for ax in axs_flat:
        ax.set_axis_off()

    #inserting each trial to a subplot
    for k in range(num_positions):
        if isinstance(axs, np.ndarray):
            ax = axs[k // cols, k % cols]
            ax.set_axis_on()

        if plot_type == 1:
            x, y = np.split(position[k],[-1], axis=1)
            ax.plot(x, y)
            ax.set_ylim([0, max(y)])
            ax.set_xlim([0, max(x)])
            if k // cols == rows-1:
                ax.set_xlabel('x')
            if k % cols == 0:
                ax.set_ylabel('y')
        elif plot_type == 2:
            ax.plot(time[k], position[k])
            if k == num_positions - 1:
                ax.legend(['x', 'y'])
            if k // cols == rows-1:
                ax.set_xlabel('time (s)')
            if k % cols == 0:
                ax.set_ylabel('position')

    plt.show()

def plot_velocity(velocity, time, plot_type=1):
    """
    Plot the velocity data against time.

    Parameters:
        velocity (list): List of velocity data arrays.
        time (list): List of time data arrays.
        plot_type (int, optional): Type of plot to generate. Default is 1.
            - plot_type = 1: Time vs. v_x & v_y velocity.
            - plot_type = 2: Time vs. tangential velocity.

    Raises:
        ValueError: If the plot_type is unknown.

    """

    if plot_type not in [1, 2]:
        raise ValueError('Unknown plot type')

    num_velocity = len(velocity)
    cols = int(np.ceil(np.sqrt(num_velocity)))
    rows = int(np.ceil(num_velocity / cols))
    #create demo figure with num subplots correlated with num of trials in data

    _ , axs = plt.subplots(rows, cols, figsize=(15, 15))
    axs_flat = axs.flatten()
    for ax in axs_flat:
        ax.set_axis_off()

    #inserting each trial to a subplot
    for k in range(num_velocity):
        if isinstance(axs, np.ndarray):
            ax = axs[k // cols, k % cols]
            ax.set_axis_on()

        if plot_type == 1:
            ax.plot(time[k], velocity[k])
            if k == num_velocity - 1:
                ax.legend(['v_x', 'v_y'])
            if k // cols == rows-1:
                ax.set_xlabel('time')
            if k % cols == 0:
                ax.set_ylabel('velocity')

        elif plot_type == 2:
            tangvel = np.sqrt(np.sum(np.square(velocity[k]), axis=1))
            ax.plot(time[k], tangvel)
            if k // cols == rows-1:
                ax.set_xlabel('time')
            if k % cols == 0:
                ax.set_ylabel('velocity')

    plt.show()

def decompose_2D(time: np.ndarray,vel: np.ndarray,
                 n_sub_movement: int = 4, x_rng: int = (-5., 5.),
                 y_rng: int = (0.1, 5)) -> tuple[float,np.ndarray,np.ndarray]:
    """
    decompose_2D - decompose two dimensional movement into submovements using the velocity profiles

    best_error, final_params, best_velocity = decompose(time,vel,numsubmovements,xrng,yrng)

    vel should be a 2 x N matrix, with the x and y velocities

    t should be a 1 x N matrix with the corresponding time (in seconds)

    n_sub_movement is the number of submovements to look for, if it is
    empty or not specified, the function will try 1 to 4 submovements

    x_rng is the valid range for the amplitude of x values (default = (-5 5))

    y_rng is the valid range for the amplitude of y values (default = (0.1 5))

    min(t0) = 0.167 * submovement number


    best_error the best (lowest) value of the error function

    best_parameters contains the function parameters corresponding to the best values
    [t0, D, Ax, Ay]. 
    If there are multiple submovements, each submovement is in different row.

    best_velocity is the velocity profile coresponding to the best values (UNIMPLANTED!!!)
    """
    # Input validation
    if time.ndim > 1:
        raise ValueError('time must be a 1D')

    if vel.ndim != 2:
        raise ValueError('vel must be an 2D ndarray - it has ' + str(vel.ndim) + ' dimensions')
    
    if vel.shape[1] != 2:
        raise ValueError('vel must be an N*2 ndarray - it is ' + str(vel.shape[0]) + '*' + str(vel.shape[1]))

    if vel.shape[0] != time.size:
        raise ValueError('vel must match time')

    # calculate tangential velocity
    tang_vel = np.sqrt(vel[:,0]**2 + vel[:,1]**2)

    lower_bounds = np.array([0,                          0.167  , x_rng[0], y_rng[0]])
    upper_bounds = np.array([max(time[-1]-0.167,0.1),    1.     , x_rng[1], y_rng[1]])
    #submovement:            start (t0)                  duration (D), Ax,  Ay

    if np.any(lower_bounds > upper_bounds):
        raise ValueError('Lower bounds exceed upper bound - infeasible')

    # initiate matrices for parameters and bounds for each submovement
    parm_per_sub = 4 # hard coded - can be change if different type of submovement is used
    init_parm = np.empty(shape=(n_sub_movement,parm_per_sub),dtype=float) # submovement parameters
    all_lower_bounds = np.empty(shape=(n_sub_movement,parm_per_sub),dtype=float) # lower bound for each parameter
    all_upper_bounds = np.empty(shape=(n_sub_movement,parm_per_sub),dtype=float) # upper bound for each parameter

    # initiate best error found
    best_error = np.inf

    # try optimazation 20 times, select the time with least error
    for _ in range(20):
        # create initial parameters for each submovement
        for iSub in range(n_sub_movement):
            init_parm[iSub,:] = lower_bounds + (upper_bounds - lower_bounds)*np.random.rand(1,parm_per_sub)
            all_upper_bounds[iSub,:] = upper_bounds.copy()
            all_lower_bounds[iSub,:] = lower_bounds.copy()
            all_lower_bounds[iSub,0] = (iSub) * 0.167

        # function to minimize
        def error_fun(params):
            epsilon = _calculate_error_MJ2D(params, time, vel, tang_vel)
            return epsilon
    
        def calculate_jacobian(params):
            jacobian = _calculate_Jacobian_MJ2D(params, time, vel, tang_vel)
            return jacobian 
    
        def calculate_hessian(params):
            hessian = _calculate_Hessian_MJ2D(params, time, vel, tang_vel)
            return hessian

        # run the optimizer - using the Hessian slows it down so don't use for now
        # Using the Jacobian speeds up the process so use it
        res = minimize(error_fun,
                       jac=calculate_jacobian,
                       #hess=calculate_hessian,
                       x0=init_parm.flatten(),
                       method='trust-constr',
                       bounds=tuple(zip(all_lower_bounds.flatten(),all_upper_bounds.flatten())),
                       options = {'maxiter':5000})

        # calculate error for the result found
        epsilon = error_fun(res.x)

        # save result if error is smaller than best found
        if epsilon < best_error:
            best_error = epsilon
            best_parameters = res.x

    # organize parameters so that every submovement is a row
    final_params = best_parameters.reshape((n_sub_movement,parm_per_sub))

    # sort by the first column (start time)
    final_params = final_params[final_params[:,0].argsort()]

    return best_error, final_params

def plot_submovements_2D(parameters, t: np.ndarray = None, plot_type: int = 1) -> tuple[plt.axes, plt.figure]:
    """
    plot_submovements_2D - plot 2D submovements after decomposition

    plot_submovements_2D(parameters,t,plot_type,x0,y0)

    The parameters should in sets of 4 for each submovement:
    [t0, D, Ax, Ay]

    plot_type:
    1 = time vs submovement velocity + sum velocity (default)
    2 = time vs submovement velocity 
    """
    if int(len(parameters))%4 != 0:
        raise ValueError('The parameters vector must have a length that is a multiple of 4')

    # parse inputs
    numsubmovements = parameters.shape[0] # each submovement is in a different row
    t0      = parameters[:, 0]
    D = parameters[:, 1]
    Ax   = parameters[:, 2]
    Ay   = parameters[:, 3]

    # make sure parameters are ordered by movement start time
    order = np.argsort(t0)
    t0      = t0[order]
    D  = D[order]
    Ax   = Ax[order]
    Ay   = Ay[order]

    # if no time was given, plot from start of first movement to end of last movement
    if t is None:
        movement_end = t0 + D # end time of each movement
        t = np.linspace(min(t0),max(movement_end),num=100)

    # init velocities
    vel_x = np.zeros((numsubmovements,t.size))
    vel_y = np.zeros((numsubmovements,t.size))

    # using minimum jerk, find velocities curve for each submovement
    for isub in range(numsubmovements):
        vel_x[isub,:], vel_y[isub,:], _ = _minimum_jerk_velocity_2D(t0[isub],D[isub], Ax [isub], Ay[isub],t)

    # get total velocity expected from submovements
    sum_vx = np.sum(vel_x,axis=0)
    sum_vy = np.sum(vel_y,axis=0)

    # create the figure
    fig, axs = plt.subplots(1, 2, figsize=(15, 5))
    ax1 = axs[0]
    ax2 = axs[1]
    vx_lines    = ax1.plot(t,vel_x.transpose(), 'b',   label=r'$V_{x}$')
    vy_lines    = ax2.plot(t,vel_y.transpose(), 'r',   label=r'$V_{y}$')
    if plot_type == 1:
        vx_sum_line = ax1.plot(t,sum_vx        , 'b--', label=r'$Sum V_{x}$')
        vy_sum_line = ax2.plot(t,sum_vy        , 'r--', label=r'$Sum V_{y}$')
        ax1.legend(handles=[vx_lines[0], vx_sum_line[0]])
        ax2.legend(handles=[vy_lines[0], vy_sum_line[0]])
    else:
        axs.legend(handles=[vx_lines[0], vy_lines[0]])

    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('x velocity')
    ax2.set_xlabel('time (s)')
    ax2.set_ylabel('y velocity')

    # return axe & figure for later plotting
    return axs, fig, vx_lines, vy_lines, vx_sum_line, vy_sum_line

def _extend_if_necessary(parameters: np.ndarray, time: np.ndarray, vel: np.ndarray, tangvel: np.ndarray, timedelta: float = 0.005) -> float:
    """
    Extends the time, velocity, and tangential velocity arrays if the last time point in the trajectory is less than the last time point in the given time array.
    Parameters:
        parameters (array): List of parameters for each submovement. Each submovement requires 4 parameters: T0, D, Dx, Dy.
        time (array): Array of time values.
        vel (array): Array of velocity values.
        tangvel (array): Array of tangential velocity values.
        timedelta (float, optional): Time interval between consecutive points. Default is 0.005.

    Returns:
        time (array): Array of time values (extended if necessary).
        vel (array): Array of velocity values (extended if necessary).
        tangvel (array): Array of tangential velocity values (extended if necessary).
    """
    # Calculate the number of submovements
    n_sub_movement = int(len(parameters)/4)

    # Find the last time point in the trajectory
    last_time = 0

    for i in range(n_sub_movement):
        t0 = parameters[i*4]
        D = parameters[i*4+1]
        last_time = max([last_time, t0+D])
    # Adjust the last time point to align with the given time interval
    last_time = (last_time*(1/timedelta))/(1/timedelta)
    # If the last time is greater than the last time point in the given time array,
    # extend the time, velocity, and tangential velocity arrays with zeros

    if last_time > time[-1]:
        new_time = np.arange(time[-1], last_time + timedelta, timedelta)
        time = np.concatenate((time[:-1], new_time))
        vel = np.concatenate((vel, np.zeros((len(time) - len(vel), vel.shape[1]))))
        tangvel = np.concatenate((tangvel, np.zeros((len(time) - len(tangvel)))))
    
    return time, vel, tangvel

def _calculate_error_MJ2D(parameters: np.ndarray,time: np.ndarray,
                         vel: np.ndarray,tangvel: np.ndarray,
                         timedelta: float = 0.005) -> float:
    """
    Calculates the error between predicted and actual trajectories in a 
    2D space based on the given parameters.
    Parameters:
        parameters (array): List of parameters for each submovement. Each submovement requires 4 parameters: T0, D, Dx, Dy.
        time (array): Array of time values.
        vel (array): Array of velocity values.
        tangvel (array): Array of tangential velocity values.
        timedelta (float, optional): Time interval between consecutive points. Default is 0.005.

    Returns:
        epsilon (float): Error between predicted and actual trajectories.
    """

    # Calculate the number of submovements
    n_sub_movement = int(len(parameters)/4)

    time, vel, tangvel = _extend_if_necessary(parameters, time, vel, tangvel, timedelta)

    # Initialize arrays for predicted trajectories
    trajectory_x = vel[:,0]
    trajectory_y = vel[:,1]

    predicted_x = np.zeros([n_sub_movement, len(time)])
    predicted_y = np.zeros([n_sub_movement, len(time)])
    predicted   = np.zeros([n_sub_movement, len(time)])

    # Calculate predicted trajectories for each submovement

    for i in range(n_sub_movement):
        t0      = parameters[i*4]
        D = parameters[i*4+1]
        Ax   = parameters[i*4+2]
        Ay   = parameters[i*4+3]

        this_rgn = np.where((time > t0) & (time < t0+D))[0]

        predicted_x[i,this_rgn], predicted_y[i,this_rgn], predicted[i,this_rgn] = \
            _minimum_jerk_velocity_2D(t0,D,Ax,Ay,time[this_rgn])

    # Calculate the sum of predicted trajectories and actual trajectories squared
    sum_predicted_x = sum(predicted_x,1)
    sum_predicted_y = sum(predicted_y,1)
    sum_predicted  = sum(predicted,1)
    sum_traj_sq = sum(trajectory_x**2 + trajectory_y**2 + tangvel**2)

    # Calculate the error between predicted and actual trajectories
    epsilon = np.sum((sum_predicted_x - trajectory_x)**2 + (sum_predicted_y - trajectory_y)**2 + (sum_predicted - tangvel)**2) / np.sum(sum_traj_sq)

    return epsilon
        
def _calculate_Jacobian_MJ2D(parameters: np.ndarray,time: np.ndarray,
                         vel: np.ndarray,tangvel: np.ndarray,
                         timedelta: float = 0.005) -> float:
    """
    Calculates the Jacobian in a 2D space based on the given parameters.
    Parameters:
        parameters (array): List of parameters for each submovement. Each submovement requires 4 parameters: T0, D, Dx, Dy.
        time (array): Array of time values.
        vel (array): Array of velocity values.
        tangvel (array): Array of tangential velocity values.
        timedelta (float, optional): Time interval between consecutive points. Default is 0.005.

    Returns:
        grad (float): Jacobian at this point 
    """
    # Calculate the number of submovements
    n_sub_movement = int(len(parameters)/4)

    trajectory_x = vel[:,0]
    trajectory_y = vel[:,1]

    predicted_x = np.zeros([n_sub_movement, len(time)])
    predicted_y = np.zeros([n_sub_movement, len(time)])
    predicted   = np.zeros([n_sub_movement, len(time)])

    J_x = np.zeros([n_sub_movement, 4*n_sub_movement, len(time)])
    J_y = np.zeros([n_sub_movement, 4*n_sub_movement, len(time)])
    J   = np.zeros([n_sub_movement, 4*n_sub_movement, len(time)])

    grad = np.zeros(4*n_sub_movement)

    for i in range(n_sub_movement):
        t0 = parameters[i*4]
        D  = parameters[i*4+1]
        Ax = parameters[i*4+2]
        Ay = parameters[i*4+3]

        this_rgn = np.where((time > t0) & (time < t0+D))[0]

        J_x[i,i*4:i*4+4,this_rgn], \
        J_y[i,i*4:i*4+4,this_rgn], \
        J[i,i*4:i*4+4,this_rgn] = \
            _minimum_jerk_Jacobian_2D(t0,D,Ax,Ay,time[this_rgn])

        predicted_x[i,this_rgn], predicted_y[i,this_rgn], predicted[i,this_rgn] = \
            _minimum_jerk_velocity_2D(t0,D,Ax,Ay,time[this_rgn])

    # Calculate the sum of predicted trajectories and actual trajectories squared
    sum_predicted_x = np.sum(predicted_x, axis=0)
    sum_predicted_y = np.sum(predicted_y, axis=0)
    sum_predicted  = np.sum(predicted, axis = 0)
    sum_traj_sq = np.sum(trajectory_x**2 + trajectory_y**2 + tangvel**2)

    sumJx = np.squeeze(np.sum(J_x, axis=0))
    sumJy = np.squeeze(np.sum(J_y, axis=0))
    sumJ = np.squeeze(np.sum(J, axis=0))
    for i in range(sumJx.shape[0]):  
        grad[i] = 2/sum_traj_sq * sum( \
            (sum_predicted_x - trajectory_x) * sumJx[i,:].T + \
            (sum_predicted_y - trajectory_y) * sumJy[i,:].T + \
            (sum_predicted - tangvel)      * sumJ[i,:].T)

    return grad


def _calculate_Hessian_MJ2D(parameters: np.ndarray,time: np.ndarray,
                         vel: np.ndarray,tangvel: np.ndarray,
                         timedelta: float = 0.005) -> float:
    """
    Calculates the Hessian in a 2D space based on the given parameters.
    Parameters:
        parameters (array): List of parameters for each submovement. Each submovement requires 4 parameters: T0, D, Dx, Dy.
        time (array): Array of time values.
        vel (array): Array of velocity values.
        tangvel (array): Array of tangential velocity values.
        timedelta (float, optional): Time interval between consecutive points. Default is 0.005.

    Returns:
        hess (float): Hessian at this point 
    """
    # Calculate the number of submovements
    n_sub_movement = int(len(parameters)/4)

    trajectory_x = vel[:,0]
    trajectory_y = vel[:,1]

    predicted_x = np.zeros([n_sub_movement, len(time)])
    predicted_y = np.zeros([n_sub_movement, len(time)])
    predicted   = np.zeros([n_sub_movement, len(time)])

    J_x = np.zeros([n_sub_movement, 4*n_sub_movement, len(time)])
    J_y = np.zeros([n_sub_movement, 4*n_sub_movement, len(time)])
    J   = np.zeros([n_sub_movement, 4*n_sub_movement, len(time)])

    H_x = np.zeros([n_sub_movement, 4*n_sub_movement, 4*n_sub_movement, len(time)])
    H_y = np.zeros([n_sub_movement, 4*n_sub_movement, 4*n_sub_movement, len(time)])
    H   = np.zeros([n_sub_movement, 4*n_sub_movement, 4*n_sub_movement, len(time)])
    
    hess = np.zeros([4*n_sub_movement, 4*n_sub_movement])

    for i in range(n_sub_movement):
        t0     =  parameters[i*4]
        D = parameters[i*4+1]
        Ax  =  parameters[i*4+2]
        Ay  =  parameters[i*4+3]

        this_rgn = np.where((time > t0) & (time < t0+D))[0]

        H_x[i,i*4:i*4+4,i*4:i*4+4,this_rgn], \
        H_y[i,i*4:i*4+4,i*4:i*4+4,this_rgn], \
        H[i,i*4:i*4+4,i*4:i*4+4, this_rgn] = \
            _minimum_jerk_Hessian_2D(t0,D,Ax,Ay,time[this_rgn])
        
        J_x[i,i*4:i*4+4,this_rgn], \
        J_y[i,i*4:i*4+4,this_rgn], \
        J[i,i*4:i*4+4,this_rgn] = \
            _minimum_jerk_Jacobian_2D(t0,D,Ax,Ay,time[this_rgn])

        predicted_x[i,this_rgn], predicted_y[i,this_rgn], predicted[i,this_rgn] = \
            _minimum_jerk_velocity_2D(t0,D,Ax,Ay,time[this_rgn])
    
    sum_predicted_x = np.sum(predicted_x, axis=0)
    sum_predicted_y = np.sum(predicted_y, axis=0)
    sum_predicted  = np.sum(predicted, axis = 0)
    sum_traj_sq = np.sum(trajectory_x**2 + trajectory_y**2 + tangvel**2)

    sumJx = np.squeeze(np.sum(J_x, axis=0))
    sumJy = np.squeeze(np.sum(J_y, axis=0))
    sumJ = np.squeeze(np.sum(J, axis=0))
    
    sumHx = np.squeeze(np.sum(H_x, axis=0))
    sumHy = np.squeeze(np.sum(H_y, axis=0))
    sumH = np.squeeze(np.sum(H, axis=0))

    for i in range(sumHx.shape[0]):
        for j in range(sumHx.shape[1]):
                hess[i,j] = 2/sum_traj_sq * sum( \
                        sumJx[i,:] * sumJx[j,:] + (sum_predicted_x - trajectory_x) * sumHx[i,j,:].T + \
                        sumJy[i,:] * sumJx[j,:] + (sum_predicted_y - trajectory_y) * sumHy[i,j,:].T + \
                        sumJ[i,:]  * sumJ[j,:]  + (sum_predicted - tangvel) * sumH[i,j,:].T)
    
    return hess

def _minimum_jerk_velocity_2D(t0: float,D: float,
                              Ax: float,Ay: float,
                              t: np.ndarray) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """
    minimumJerkVelocity2D - evaluate a minimum jerk velocity curve with separate displacement for x / y

    see Flash and Hogan (1985) for details on the minimum jerk equation

        t0 = movement start time (scalar)
        D  = movement duration (scalar)
        Ax   = displacement resulting from the movement (x) (scalar)
        Ay   = displacement resulting from the movement (y) (scalar)

    The function is evaluated at times t (vector)

    x_vel, y_vel and tan_vel are the x velocity, y velocity and tangential velocities
    """
    # normalise time to t0 and movement duration, take only the time of the movement
    normlized_time = (t - t0)/D
    logical_movement = (normlized_time >= 0) & (normlized_time <= 1)

    # normalise displacement to movement duration
    norm_disp_x = Ax/D
    norm_disp_y = Ay/D
    tang_norm_disp = np.sqrt(norm_disp_x**2 + norm_disp_y**2)

    # make x_vel, y_vel & tan_vel that are zero outside of calculated area
    x_vel = np.zeros(t.size)
    y_vel = np.zeros(t.size)
    tan_vel  = np.zeros(t.size)

    # calculate velocities
    def min_jerk_2d_fun(base_val):
        # the polynomial function from Flash and Hogan (1985)
        return base_val * (-60*normlized_time[logical_movement]**3 + 30*normlized_time[logical_movement]**4 + 30*normlized_time[logical_movement]**2)

    x_vel[logical_movement]   = min_jerk_2d_fun(norm_disp_x)
    y_vel[logical_movement]   = min_jerk_2d_fun(norm_disp_y)
    tan_vel[logical_movement] = min_jerk_2d_fun(tang_norm_disp)

    return x_vel, y_vel, tan_vel

def _minimum_jerk_Jacobian_2D(t0: float,D: float,
                              Ax: float,Ay: float,
                              t: np.ndarray) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """
    minimumJerkJacobian2D - evaluate the Jacobian (partial derivative) of a minimum jerk velocity curve with separate displacement for x / y

        t0 = movement start time (scalar)
        D  = movement duration (scalar)
        Ax   = displacement resulting from the movement (x) (scalar)
        Ay   = displacement resulting from the movement (y) (scalar)

    The function is evaluated at times t (vector)

    J_x, J_y and J are the gradients (partial derivatives) of the same quantities

    For the derivation, see https://github.com/JasonFriedman/submovements/tree/master/matlab_deriveGradiantHessian 
    """
    # normalise time to t0 and movement duration, take only the time of the movement
    normlized_time = (t - t0)/D
    logical_movement = (normlized_time >= 0) & (normlized_time <= 1)

    J_x = np.zeros((t.size,4))
    J_y = np.zeros((t.size,4))
    J   = np.zeros((t.size,4))

    J_x[logical_movement, 0] = -Ax * ((1.0 / D**3 * (t[logical_movement]*2.0 - t0*2.0) * 30) - (1.0 / D**4 * (t[logical_movement] - t0)**2 * 180) + (1.0 / D**5 * (t[logical_movement] - t0)**3 * 120))
    J_x[logical_movement, 1] = -Ax * ((1.0 / D**4 * (t[logical_movement] - t0)**2 * 90) - (1.0 / D**5 * (t[logical_movement] - t0)**3 * 240) + (1.0 / D**6 * (t[logical_movement] - t0)**4 * 150))
    J_x[logical_movement, 2] = (1.0 / D**3 * (t[logical_movement] - t0)**2 * 30) - (1.0 / D**4 * (t[logical_movement] - t0)**3 * 60) + (1.0 / D**5 * (t[logical_movement] - t0)**4 * 30)

    J_y[logical_movement, 0] = -Ay * ((1.0 / D**3 * (t[logical_movement]*2.0 - t0*2.0) * 30) - (1.0 / D**4 * (t[logical_movement] - t0)**2 * 180) + (1.0 / D**5 * (t[logical_movement] - t0)**3 * 120))
    J_y[logical_movement, 1] = -Ay * ((1.0 / D**4 * (t[logical_movement] - t0)**2 * 90) - (1.0 / D**5 * (t[logical_movement] - t0)**3 * 240) + (1.0 / D**6 * (t[logical_movement] - t0)**4 * 150))
    J_y[logical_movement, 2] = (1.0 / D**3 * (t[logical_movement] - t0)**2 * 30) - (1.0 / D**4 * (t[logical_movement] - t0)**3 * 60) + (1.0 / D**5 * (t[logical_movement] - t0)**4 * 30)

    J[logical_movement, 0] = ((1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**3 * (D - t[logical_movement] + t0)**4 * 4.0 - 1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**3 * 4.0) / np.sqrt(1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4) * -15)
    J[logical_movement, 1] = ((1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**3 * 4.0 - 1.0 / D**11 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4 * 10) / np.sqrt(1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4) * 15)
    J[logical_movement, 2] = Ax * 1.0 / D**10 * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4 / np.sqrt(1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4) * 30
    J[logical_movement, 3] = Ay * 1.0 / D**10 * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4 / np.sqrt(1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4) * 30
    
    return J_x, J_y, J

def _minimum_jerk_Hessian_2D(t0: float,D: float,
                              Ax: float,Ay: float,
                              t: np.ndarray) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """
    minimumJerkHessian2D - evaluate the Hessian (second-order partial derivatives) of a minimum jerk velocity curve with separate displacement for x / y

        t0 = movement start time (scalar)
        D  = movement duration (scalar)
        Ax   = displacement resulting from the movement (x) (scalar)
        Ay   = displacement resulting from the movement (y) (scalar)

    The function is evaluated at times t (vector)

    H_x, H_y and H are the Hessian (second-order partial derivatives) 

    For the derivation, see https://github.com/JasonFriedman/submovements/tree/master/matlab_deriveGradiantHessian
    """
    # normalise time to t0 and movement duration, take only the time of the movement
    normlized_time = (t - t0)/D
    logical_movement = (normlized_time >= 0) & (normlized_time <= 1)

    H_x = np.zeros((t.size,4,4))
    H_y = np.zeros((t.size,4,4))
    H   = np.zeros((t.size,4,4))

    H_x[logical_movement, 0, 0] = Ax * ((1.0 / D**4 * (t[logical_movement]*2.0 - t0*2.0) * -180) + (1.0 / D**5 * (t[logical_movement] - t0)**2 * 360) + (1.0 / D**3 * 60))
    H_x[logical_movement, 0, 1] = Ax * ((1.0 / D**4 * (t[logical_movement]*2.0 - t0*2.0) * 90) - (1.0 / D**5 * (t[logical_movement] - t0)**2 * 720) + (1.0 / D**6 * (t[logical_movement] - t0)**3 * 600))
    H_x[logical_movement, 0, 2] = (1.0 / D**3 * (t[logical_movement]*2.0 - t0*2.0) * -30) + (1.0 / D**4 * (t[logical_movement] - t0)**2 * 180) - (1.0 / D**5 * (t[logical_movement] - t0)**3 * 120)

    H_x[logical_movement, 1, 0] = Ax * ((1.0 / D**4 * (t[logical_movement]*2.0 - t0*2.0) * 90) - (1.0 / D**5 * (t[logical_movement] - t0)**2 * 720) + (1.0 / D**6 * (t[logical_movement] - t0)**3 * 600))
    H_x[logical_movement, 1, 1] = Ax * ((1.0 / D**5 * (t[logical_movement] - t0)**2 * 360) - (1.0 / D**6 * (t[logical_movement] - t0)**3 * 1200) + (1.0 / D**7 * (t[logical_movement] - t0)**4 * 900))
    H_x[logical_movement, 1, 2] = (1.0 / D**4 * (t[logical_movement] - t0)**2 * -90) + (1.0 / D**5 * (t[logical_movement] - t0)**3 * 240) - (1.0 / D**6 * (t[logical_movement] - t0)**4 * 150)

    H_x[logical_movement, 2, 0] = (1.0 / D**3 * (t[logical_movement]*2.0 - t0*2.0) * -30) + (1.0 / D**4 * (t[logical_movement] - t0)**2 * 180) - (1.0 / D**5 * (t[logical_movement] - t0)**3 * 120)
    H_x[logical_movement, 2, 1] = (1.0 / D**4 * (t[logical_movement] - t0)**2 * -90) + (1.0 / D**5 * (t[logical_movement] - t0)**3 * 240) - (1.0 / D**6 * (t[logical_movement] - t0)**4 * 150)

    H_y[logical_movement, 0, 0] = Ay * ((1.0 / D**4 * (t[logical_movement]*2.0 - t0*2.0) * -180) + (1.0 / D**5 * (t[logical_movement] - t0)**2 * 360) + (1.0 / D**3 * 60))
    H_y[logical_movement, 0, 1] = Ay * ((1.0 / D**4 * (t[logical_movement]*2.0 - t0*2.0) * 90) - (1.0 / D**5 * (t[logical_movement] - t0)**2 * 720) + (1.0 / D**6 * (t[logical_movement] - t0)**3 * 600))

    H_y[logical_movement, 0, 3] = (1.0 / D**3 * (t[logical_movement]*2.0 - t0*2.0) * -30) + (1.0 / D**4 * (t[logical_movement] - t0)**2 * 180) - (1.0 / D**5 * (t[logical_movement] - t0)**3 * 120)
    H_y[logical_movement, 1, 0] = Ay * ((1.0 / D**4 * (t[logical_movement]*2.0 - t0*2.0) * 90) - (1.0 / D**5 * (t[logical_movement] - t0)**2 * 720) + (1.0 / D**6 * (t[logical_movement] - t0)**3 * 600))
    H_y[logical_movement, 1, 1] = Ay * ((1.0 / D**5 * (t[logical_movement] - t0)**2 * 360) - (1.0 / D**6 * (t[logical_movement] - t0)**3 * 1200) + (1.0 / D**7 * (t[logical_movement] - t0)**4 * 900))

    H_y[logical_movement, 1, 3] = (1.0 / D**4 * (t[logical_movement] - t0)**2 * -90) + (1.0 / D**5 * (t[logical_movement] - t0)**3 * 240) - (1.0 / D**6 * (t[logical_movement] - t0)**4 * 150)

    H_y[logical_movement, 3, 0] = (1.0 / D**3 * (t[logical_movement]*2.0 - t0*2.0) * -30) + (1.0 / D**4 * (t[logical_movement] - t0)**2 * 180) - (1.0 / D**5 * (t[logical_movement] - t0)**3 * 120)
    H_y[logical_movement, 3, 1] = (1.0 / D**4 * (t[logical_movement] - t0)**2 * -90) + (1.0 / D**5 * (t[logical_movement] - t0)**3 * 240) - (1.0 / D**6 * (t[logical_movement] - t0)**4 * 150)

    H[logical_movement, 0, 0] = 1.0 / D**20 * (Ax**2 + Ay**2)**2 * (t[logical_movement] - t0)**6 * (D - t[logical_movement] + t0)**6 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * (D * t[logical_movement] * -6.0 + D * t0 * 6.0 - t[logical_movement] * t0 * 12.0 + D**2 + t[logical_movement]**2 * 6.0 + t0**2 * 6.0) * 60.0
    H[logical_movement, 0, 1] = 1.0 / D**21 * (Ax**2 + Ay**2)**2 * (t[logical_movement] - t0)**7 * (D - t[logical_movement] + t0)**6 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * (D * t[logical_movement] * -12.0 + D * t0 * 12.0 - t[logical_movement] * t0 * 20.0 + D**2 * 3.0 + t[logical_movement]**2 * 10.0 + t0**2 * 10.0) * 60.0
    H[logical_movement, 0, 2] = Ax / D**20 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**7 * (D - t[logical_movement] * 2.0 + t0 * 2.0) * (D - t[logical_movement] + t0)**7 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * -60.0
    H[logical_movement, 0, 3] = Ay / D**20 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**7 * (D - t[logical_movement] * 2.0 + t0 * 2.0) * (D - t[logical_movement] + t0)**7 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * -60.0
    H[logical_movement, 1, 0] = 1.0 / D**21 * (Ax**2 + Ay**2)**2 * (t[logical_movement] - t0)**7 * (D - t[logical_movement] + t0)**6 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * (D * t[logical_movement] * -12 + D * t0 * 12 - t[logical_movement] * t0 * 20 + D**2 * 3 + t[logical_movement]**2 * 10 + t0**2 * 10) * 60
    H[logical_movement, 1, 1] = 1.0 / D**22 * (Ax**2 + Ay**2)**2 * (t[logical_movement] - t0)**8 * (D - t[logical_movement] + t0)**6 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * (D * t[logical_movement] * -20 + D * t0 * 20 - t[logical_movement] * t0 * 30 + D**2 * 6 + t[logical_movement]**2 * 15 + t0**2 * 15) * 60
    H[logical_movement, 1, 2] = Ax / D**21 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**8 * (D * 3 - t[logical_movement] * 5 + t0 * 5) * (D - t[logical_movement] + t0)**7 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * -30
    H[logical_movement, 1, 3] = Ay / D**21 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**8 * (D * 3 - t[logical_movement] * 5 + t0 * 5) * (D - t[logical_movement] + t0)**7 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * -30
    H[logical_movement, 2, 0] = Ax / D**20 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**7 * (D - t[logical_movement]*2.0 + t0*2.0) * (D - t[logical_movement] + t0)**7 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * -60.0
    H[logical_movement, 2, 1] = Ax / D**21 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**8 * (D*3.0 - t[logical_movement]*5.0 + t0*5.0) * (D - t[logical_movement] + t0)**7 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * -30.0
    H[logical_movement, 2, 2] = Ay**2 / D**20 * (t[logical_movement] - t0)**8 * (D - t[logical_movement] + t0)**8 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * 30.0
    H[logical_movement, 2, 3] = Ax * Ay / D**20 * (t[logical_movement] - t0)**8 * (D - t[logical_movement] + t0)**8 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * -30.0
    H[logical_movement, 3, 0] = Ay * 1.0 / D**20 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**7 * (D - t[logical_movement]*2.0 + t0*2.0) * (D - t[logical_movement] + t0)**7 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * -60.0
    H[logical_movement, 3, 1] = Ay * 1.0 / D**21 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**8 * (D*3.0 - t[logical_movement]*5.0 + t0*5.0) * (D - t[logical_movement] + t0)**7 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * -30.0
    H[logical_movement, 3, 2] = Ax * Ay * 1.0 / D**20 * (t[logical_movement] - t0)**8 * (D - t[logical_movement] + t0)**8 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * -30.0
    H[logical_movement, 3, 3] = Ax**2 * 1.0 / D**20 * (t[logical_movement] - t0)**8 * (D - t[logical_movement] + t0)**8 / (1.0 / D**10 * (Ax**2 + Ay**2) * (t[logical_movement] - t0)**4 * (D - t[logical_movement] + t0)**4)**(3.0 / 2.0) * 30.0 
    
    return H_x, H_y, H 