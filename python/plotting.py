"""Plot helpers for movement trajectories and decomposed submovements."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from kinematics import prepare_submovement_kinematics


def plot_position(position, time, plot_type=1):
    """Plot one or more position trajectories.

    Parameters
    ----------
    position : list of ndarray
        Sequence of position arrays, one per movement. Each array should be
        shaped ``(N, 2)`` for x/y position.
    time : list of ndarray
        Sequence of time vectors aligned with ``position``.
    plot_type : int, default=1
        Plot mode:
        1. x versus y trajectory.
        2. x and y versus time.

    Returns
    -------
    None
        The function draws the figures and calls ``plt.show()``.
    """
    if plot_type not in [1, 2]:
        raise ValueError('Unknown plot type')

    n = len(position)
    cols = int(np.ceil(np.sqrt(n)))
    rows = int(np.ceil(n / cols))
    _, axs = plt.subplots(rows, cols, figsize=(15, 15))
    axs = np.atleast_1d(axs).ravel()
    for ax in axs:
        ax.set_axis_off()

    for k in range(n):
        ax = axs[k]
        ax.set_axis_on()
        if plot_type == 1:
            ax.plot(position[k][:, 0], position[k][:, 1])
            if k // cols == rows - 1:
                ax.set_xlabel('x')
            if k % cols == 0:
                ax.set_ylabel('y')
        else:
            ax.plot(time[k], position[k])
            if k == n - 1:
                ax.legend(['x', 'y'])
            if k // cols == rows - 1:
                ax.set_xlabel('time (s)')
            if k % cols == 0:
                ax.set_ylabel('position')

    plt.show()


def plot_velocity(velocity, time, plot_type=1):
    """Plot one or more velocity trajectories.

    Parameters
    ----------
    velocity : list of ndarray
        Sequence of velocity arrays, one per movement. Each array should be
        shaped ``(N, 2)`` for x/y velocity.
    time : list of ndarray
        Sequence of time vectors aligned with ``velocity``.
    plot_type : int, default=1
        Plot mode:
        1. x and y velocity versus time.
        2. Tangential velocity versus time.

    Returns
    -------
    None
        The function draws the figures and calls ``plt.show()``.
    """
    if plot_type not in [1, 2]:
        raise ValueError('Unknown plot type')

    n = len(velocity)
    cols = int(np.ceil(np.sqrt(n)))
    rows = int(np.ceil(n / cols))
    _, axs = plt.subplots(rows, cols, figsize=(15, 15))
    axs = np.atleast_1d(axs).ravel()
    for ax in axs:
        ax.set_axis_off()

    for k in range(n):
        ax = axs[k]
        ax.set_axis_on()
        if plot_type == 1:
            ax.plot(time[k], velocity[k])
            if k == n - 1:
                ax.legend(['v_x', 'v_y'])
        else:
            tang = np.linalg.norm(velocity[k], axis=1)
            ax.plot(time[k], tang)
        if k // cols == rows - 1:
            ax.set_xlabel('time')
        if k % cols == 0:
            ax.set_ylabel('velocity')

    plt.show()


def plot_submovements_1d(parameters, t=None, plot_type=1, x0=0):
    """Plot a 1D submovement decomposition.

    Parameters
    ----------
    parameters : array-like
        Flattened or 2D parameter array with one row per submovement in the
        order ``[t0, D, A]``.
    t : ndarray, optional
        Time vector used for plotting. If omitted, a default plotting grid is
        generated from the submovement parameters.
    plot_type : int, default=1
        Plot mode:
        1. Individual submovement velocities and summed velocity.
        2. Individual submovement velocities only.
        3. Individual submovement positions and summed position.
        4. Individual submovement positions only.
    x0 : float, default=0
        Initial position offset.

    Returns
    -------
    None
        The function plots onto the current figure.
    """
    data = prepare_submovement_kinematics(parameters, 1, t, [x0])
    t = data['t']
    t0 = data['t0']
    duration = data['D']
    v = data['vel'][:, :, 0]
    x = data['pos'][:, :, 0]
    x_sum = data['sumPos'][:, 0]

    if plot_type in (1, 2):
        h = plt.plot(t, v.T, 'b')
        plt.xlabel('time')
        plt.ylabel('velocity')
        if plot_type == 1:
            hh = plt.plot(t, np.sum(v, axis=0), 'k--', linewidth=2)
            plt.legend([h[0], hh[0]], ['Submovements velocity', 'Sum movements velocity'])
    elif plot_type in (3, 4):
        h = None
        for s in range(x.shape[0]):
            r = (t >= t0[s]) & (t <= t0[s] + duration[s])
            lines = plt.plot(t[r], x[s, r], 'b')
            h = lines if h is None else h
        plt.xlabel('time')
        plt.ylabel('position')
        if plot_type == 3:
            hh = plt.plot(t, x_sum, 'k--', linewidth=2)
            plt.legend([h[0], hh[0]], ['Submovements x', 'Sum submovements x'])
        else:
            plt.legend([h[0]], ['Submovements x'])
    else:
        raise ValueError('Unknown plot type')


def plot_submovements_2d(parameters, t=None, plot_type=1, x0=0, y0=0):
    """Plot a 2D submovement decomposition.

    Parameters
    ----------
    parameters : array-like
        Flattened or 2D parameter array with one row per submovement in the
        order ``[t0, D, Ax, Ay]``.
    t : ndarray, optional
        Time vector used for plotting. If omitted, a default plotting grid is
        generated from the submovement parameters.
    plot_type : int, default=1
        Plot mode:
        1. Side-by-side velocity plots for x and y, including summed velocity.
        2. Side-by-side velocity plots for x and y, submovements only.
        3. Side-by-side position plots for x and y, including summed position.
        4. Side-by-side position plots for x and y, submovements only.
        5. Spatial x-y trajectory plot.
    x0 : float, default=0
        Initial x-position offset.
    y0 : float, default=0
        Initial y-position offset.

    Returns
    -------
    numpy.ndarray | matplotlib.axes.Axes
        Returns a two-element axes array for plot types 1-4, and a single axes
        object for plot type 5.
    """
    data = prepare_submovement_kinematics(parameters, 2, t, [x0, y0])
    t = data['t']
    t0 = data['t0']
    duration = data['D']
    vx = data['vel'][:, :, 0]
    vy = data['vel'][:, :, 1]
    x = data['pos'][:, :, 0]
    y = data['pos'][:, :, 1]
    x_sum = data['sumPos'][:, 0]
    y_sum = data['sumPos'][:, 1]

    if plot_type in (1, 2):
        _, axs = plt.subplots(1, 2, figsize=(12, 4), sharex=True)

        hx = axs[0].plot(t, vx.T, 'b')
        axs[0].set_xlabel('time')
        axs[0].set_ylabel('velocity x')

        hy = axs[1].plot(t, vy.T, 'r')
        axs[1].set_xlabel('time')
        axs[1].set_ylabel('velocity y')

        if plot_type == 1:
            hxs = axs[0].plot(t, np.sum(vx, axis=0), 'k--', linewidth=2)
            hys = axs[1].plot(t, np.sum(vy, axis=0), 'k--', linewidth=2)
            axs[0].legend([hx[0], hxs[0]], ['Submovements v_x', 'Sum movements v_x'])
            axs[1].legend([hy[0], hys[0]], ['Submovements v_y', 'Sum movements v_y'])

        plt.tight_layout()
        return axs
    elif plot_type in (3, 4):
        _, axs = plt.subplots(1, 2, figsize=(12, 4), sharex=True)

        hx = None
        hy = None
        for s in range(x.shape[0]):
            r = (t >= t0[s]) & (t <= t0[s] + duration[s])
            lines_x = axs[0].plot(t[r], x[s, r], 'b')
            lines_y = axs[1].plot(t[r], y[s, r], 'r')
            hx = lines_x if hx is None else hx
            hy = lines_y if hy is None else hy

        axs[0].set_xlabel('time')
        axs[0].set_ylabel('position x')
        axs[1].set_xlabel('time')
        axs[1].set_ylabel('position y')

        if plot_type == 3:
            hxs = axs[0].plot(t, x_sum, 'k--', linewidth=2)
            hys = axs[1].plot(t, y_sum, 'k--', linewidth=2)
            axs[0].legend([hx[0], hxs[0]], ['Submovements x', 'Sum x'])
            axs[1].legend([hy[0], hys[0]], ['Submovements y', 'Sum y'])

        plt.tight_layout()
        return axs
    elif plot_type == 5:
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        h = None
        for s in range(x.shape[0]):
            r = (t >= t0[s]) & (t <= t0[s] + duration[s])
            lines = ax.plot(x[s, r], y[s, r], 'b')
            h = lines if h is None else h
        hh = ax.plot(x_sum, y_sum, 'k--', linewidth=2)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.legend([h[0], hh[0]], ['Submovements', 'Sum submovements'])
        fig.tight_layout()
        return ax
    else:
        raise ValueError('Unknown plot type')


def plot_velocity_and_reconstruction_2d(velocity, parameters, t=None, x0=0, y0=0):
    """Plot measured and reconstructed 2D velocity side by side.

    Parameters
    ----------
    velocity : ndarray
        Measured 2D velocity array shaped ``(N, 2)``.
    parameters : array-like
        Flattened or 2D parameter array with one row per submovement in the
        order ``[t0, D, Ax, Ay]``.
    t : ndarray, optional
        Time vector aligned with ``velocity``. If omitted, a default plotting
        grid is generated from the submovement parameters.
    x0 : float, default=0
        Initial x-position offset used when preparing the reconstruction.
    y0 : float, default=0
        Initial y-position offset used when preparing the reconstruction.

    Returns
    -------
    numpy.ndarray
        A two-element axes array: left axis for x velocity and right axis for
        y velocity.
    """
    velocity = np.asarray(velocity, dtype=float)
    if velocity.ndim != 2 or velocity.shape[1] != 2:
        raise ValueError('velocity must have shape (N, 2)')

    data = prepare_submovement_kinematics(parameters, 2, t, [x0, y0])
    t = data['t']
    reconstructed = np.sum(data['vel'], axis=0)

    if velocity.shape[0] != t.shape[0]:
        raise ValueError('velocity and t must have matching lengths')

    _, axs = plt.subplots(1, 2, figsize=(12, 4), sharex=True)

    axs[0].plot(t, velocity[:, 0], 'k', label='Measured v_x')
    axs[0].plot(t, reconstructed[:, 0], 'b--', linewidth=2, label='Reconstructed v_x')
    axs[0].set_xlabel('time')
    axs[0].set_ylabel('velocity x')
    axs[0].legend()

    axs[1].plot(t, velocity[:, 1], 'k', label='Measured v_y')
    axs[1].plot(t, reconstructed[:, 1], 'r--', linewidth=2, label='Reconstructed v_y')
    axs[1].set_xlabel('time')
    axs[1].set_ylabel('velocity y')
    axs[1].legend()

    plt.tight_layout()
    return axs


def plot_submovements_3d(parameters, t=None, plot_type=1, x0=0, y0=0, z0=0):
    """Plot a 3D submovement decomposition.

    Parameters
    ----------
    parameters : array-like
        Flattened or 2D parameter array with one row per submovement in the
        order ``[t0, D, Ax, Ay, Az]``.
    t : ndarray, optional
        Time vector used for plotting. If omitted, a default plotting grid is
        generated from the submovement parameters.
    plot_type : int, default=1
        Plot mode:
        1. Component velocities and summed velocities versus time.
        2. Component velocities only.
        3. Component positions and summed positions versus time.
        4. Component positions only.
        5. Spatial x-y trajectory plot of the decomposition.
    x0 : float, default=0
        Initial x-position offset.
    y0 : float, default=0
        Initial y-position offset.
    z0 : float, default=0
        Initial z-position offset.

    Returns
    -------
    None
        The function plots onto the current figure.
    """
    data = prepare_submovement_kinematics(parameters, 3, t, [x0, y0, z0])
    t = data['t']
    t0 = data['t0']
    duration = data['D']
    vx = data['vel'][:, :, 0]
    vy = data['vel'][:, :, 1]
    vz = data['vel'][:, :, 2]
    x = data['pos'][:, :, 0]
    y = data['pos'][:, :, 1]
    z = data['pos'][:, :, 2]
    x_sum = data['sumPos'][:, 0]
    y_sum = data['sumPos'][:, 1]
    z_sum = data['sumPos'][:, 2]

    if plot_type in (1, 2):
        h1 = plt.plot(t, vx.T, 'b')
        h2 = plt.plot(t, vy.T, 'r')
        h3 = plt.plot(t, vz.T, 'g')
        plt.xlabel('time')
        plt.ylabel('velocity')
        if plot_type == 1:
            hh1 = plt.plot(t, np.sum(vx, axis=0), 'k--', linewidth=2)
            hh2 = plt.plot(t, np.sum(vy, axis=0), 'm--', linewidth=2)
            hh3 = plt.plot(t, np.sum(vz, axis=0), 'c--', linewidth=2)
            plt.legend([h1[0], h2[0], h3[0], hh1[0], hh2[0], hh3[0]], ['Submovements v_x', 'Submovements v_y', 'Submovements v_z', 'Sum v_x', 'Sum v_y', 'Sum v_z'])
    elif plot_type in (3, 4):
        h = None
        for s in range(x.shape[0]):
            r = (t >= t0[s]) & (t <= t0[s] + duration[s])
            lines = plt.plot(t[r], x[s, r], 'b')
            plt.plot(t[r], y[s, r], 'r')
            plt.plot(t[r], z[s, r], 'g')
            h = lines if h is None else h
        plt.xlabel('time')
        plt.ylabel('position')
        if plot_type == 3:
            hh1 = plt.plot(t, x_sum, 'k--', linewidth=2)
            hh2 = plt.plot(t, y_sum, 'm--', linewidth=2)
            hh3 = plt.plot(t, z_sum, 'c--', linewidth=2)
            plt.legend([h[0], hh1[0], hh2[0], hh3[0]], ['Submovements x', 'Sum x', 'Sum y', 'Sum z'])
    elif plot_type == 5:
        h = None
        for s in range(x.shape[0]):
            r = (t >= t0[s]) & (t <= t0[s] + duration[s])
            lines = plt.plot(x[s, r], y[s, r], 'b')
            h = lines if h is None else h
        hh = plt.plot(x_sum, y_sum, 'k--', linewidth=2)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend([h[0], hh[0]], ['Submovements', 'Sum submovements'])
    else:
        raise ValueError('Unknown plot type')
