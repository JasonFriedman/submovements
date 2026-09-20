# submovements
This is code for the decomposition of velocity data into submovements

There are two versions: one for [Matlab](matlab) and one for [Python](python)

## Matlab
The matlab code requires the Optimization toolkit to be installed
The file [sample.m](matlab/sample.m) gives an example of how to use the toolkit

For long 2D trajectories, you can fit in windows using `decompose2Dwindows`:

```matlab
time = (0:0.005:20)';
vel = yourVelocityMatrix; % N x 2

[bestErrors,bestParameters,bestVelocity,decomposition] = decompose2Dwindows(...
	time,vel,1:4,[-5 5],[0.1 5],0.03,3);

% decomposition fields include:
% t0s, Ds, Axs, Ays, endtimes, submovementsVelocity, reconstructedVelocity
```

### 1D decomposition examples
`decompose1D` fits one-dimensional velocity data with minimum-jerk submovements.

```matlab
time = (0:0.005:1.5)';
vel = yourVelocityVector; % N x 1

% Try 1 to 4 submovements (default range), using default amplitude bounds
[bestError,bestParameters,bestVelocity] = decompose1D(time,vel);

% Explicitly fit a fixed number of submovements with custom amplitude range
[bestError2,bestParameters2,bestVelocity2] = decompose1D(time,vel,2,[-3 3]);

% Parameters are ordered as [t0 D A] per submovement
```

`decompose1Dwindows` is useful for long trials by fitting consecutive time windows.

```matlab
time = (0:0.005:20)';
vel = yourVelocityVector; % N x 1

submovementRange = 1:4;
criteria = 0.03;
windowSize = 3; % seconds

[bestErrors,bestParameters,bestVelocity,decomposition] = decompose1Dwindows(...
	time,vel,submovementRange,[-5 5],criteria,windowSize);

% decomposition fields include:
% t0s, Ds, As, endtimes, submovementsVelocity, reconstructedVelocity
```

## Python
The [Jupyter](python/run_example.ipynb) notebook gives an end-to-end example.

Main Python modules:

- `decompose.py` — `decompose_1d`, `decompose_2d`, `decompose_3d`
- `windows.py` — adaptive windowed decomposition (`decompose_*_windows`)
- `plotting.py` — position/velocity and submovement plotting helpers
- `metrics.py` — overlap and relative-onset metrics
- `experimental_data_loader.py` — experiment-specific CSV loading utilities

The compatibility module `movement_decompose_2d.py` still exposes legacy names (`decompose_2D`, `plot_submovements_2D`, etc.) while routing to the new modular implementation.

## Background
An older version of the (Matlab) code with more explanation can be found [here](https://noisyaccumulation.blogspot.com/2012/02/how-to-decompose-2d-trajectory-data.html)
