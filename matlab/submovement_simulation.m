function varargout = submovement_simulation(varargin)
% SUBMOVEMENT_SIMULATION MATLAB code for submovement_simulation.fig
%      SUBMOVEMENT_SIMULATION, by itself, creates a new SUBMOVEMENT_SIMULATION or raises the existing
%      singleton*.
%
%      H = SUBMOVEMENT_SIMULATION returns the handle to a new SUBMOVEMENT_SIMULATION or the handle to
%      the existing singleton*.
%
%      SUBMOVEMENT_SIMULATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SUBMOVEMENT_SIMULATION.M with the given input arguments.
%
%      SUBMOVEMENT_SIMULATION('Property','Value',...) creates a new SUBMOVEMENT_SIMULATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before submovement_simulation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to submovement_simulation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help submovement_simulation

% Last Modified by GUIDE v2.5 09-Jun-2021 23:57:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @submovement_simulation_OpeningFcn, ...
                   'gui_OutputFcn',  @submovement_simulation_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before submovement_simulation is made visible.
function submovement_simulation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to submovement_simulation (see VARARGIN)

% Choose default command line output for submovement_simulation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

plotgraph(handles)

% UIWAIT makes submovement_simulation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = submovement_simulation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function plotgraph(handles)

axes(handles.axes1)
cla

A = str2double(handles.A.String);
D = str2double(handles.D.String);
t0 = str2double(handles.t0.String);
sigma = str2double(handles.sigma.String);
mu = str2double(handles.mu.String);

A2 = str2double(handles.A2.String);
D2 = str2double(handles.D2.String);
t02 = str2double(handles.t02.String);
sigma2 = str2double(handles.sigma2.String);
mu2 = str2double(handles.mu2.String);


if handles.onesubmovement.Value
    handles.A2Label.Visible = 'off';
    handles.A2.Visible = 'off';
    handles.D2Label.Visible = 'off';
    handles.D2.Visible = 'off';
    handles.t02Label.Visible = 'off';
    handles.t02.Visible = 'off';
    handles.sigma2Label.Visible = 'off';
    handles.sigma2.Visible = 'off'; 
    handles.mu2Label.Visible = 'off';
    handles.mu2.Visible = 'off';
    handles.x02Label.Visible = 'off';
    handles.x02.Visible = 'off';
else
    handles.A2Label.Visible = 'on';
    handles.A2.Visible = 'on';
    handles.D2Label.Visible = 'on';
    handles.D2.Visible = 'on';
    handles.t02Label.Visible = 'on';
    handles.t02.Visible = 'on';
end


if isnan(A) || isnan(D) || isnan(t0)
    return
end

if handles.position.Value
    handles.x0.Visible = 'on';
    handles.x0Label.Visible = 'on';
    if handles.twosubmovements.Value
        handles.x02.Visible = 'on';
        handles.x02Label.Visible = 'on';
    end

    x0 = str2double(handles.x0.String);
    x02 = str2double(handles.x02.String);
    if isnan(x0)
        return
    end
else
    handles.x0.Visible = 'off';
    handles.x0Label.Visible = 'off';
    handles.x02.Visible = 'off';
    handles.x02Label.Visible = 'off';
end

if handles.onesubmovement.Value
    t = linspace(t0,t0+D,100);
else
    t = linspace(min([t0 t02]),max([t0+D t02+D2]));
end

if handles.minimumJerk.Value
    handles.mu.Visible = 'off';
    handles.muLabel.Visible = 'off';
    handles.sigma.Visible = 'off';
    handles.sigmaLabel.Visible = 'off';
    handles.sigma2.Visible = 'off';
    handles.sigma2Label.Visible = 'off';
    handles.mu2.Visible = 'off';
    handles.mu2Label.Visible = 'off';

    if handles.velocity.Value
        y = minimumJerkVelocity1D(t0,D,A,t);
    else
        y = minimumJerkPosition1D(t0,D,A,x0,t);
    end
    if handles.twosubmovements.Value
        if handles.velocity.Value
            y2 = minimumJerkVelocity1D(t02,D2,A2,t);
        else
            y2 = minimumJerkPosition1D(t02,D2,A2,x02,t);
        end
    end

else
    handles.mu.Visible = 'on';
    handles.muLabel.Visible = 'on';
    handles.sigma.Visible = 'on';
    handles.sigmaLabel.Visible = 'on';
    if handles.twosubmovements.Value
        handles.mu2.Visible = 'on';
        handles.mu2Label.Visible = 'on';
        handles.sigma2.Visible = 'on';
        handles.sigma2Label.Visible = 'on';
    end

    if handles.velocity.Value
        y = lgnbVelocity1D(t0,D,A,sigma,mu,t);
    else
        y = lgnbPosition1D(t0,D,A,sigma,mu,x0,t);
    end
    if handles.twosubmovements.Value
        if handles.velocity.Value
            y2 = lgnbVelocity1D(t02,D2,A2,sigma2,mu2,t);
        else
            y2 = lgnbPosition1D(t02,D2,A2,sigma2,mu2,x02,t);
        end
    end

end

plot(t,y,'LineWidth',2);
if handles.twosubmovements.Value
    hold on;
    plot(t,y2,'LineWidth',2);
    plot(t,y+y2,'--','LineWidth',2);
end
xlabel('time');
if handles.position.Value
    ylabel('Position');
else
    ylabel('Velocity');
end
set(gca,'FontSize',14)

function t0_Callback(hObject, eventdata, handles)
plotgraph(handles);

function D_Callback(hObject, eventdata, handles)
plotgraph(handles);

function A_Callback(hObject, eventdata, handles)
plotgraph(handles);

function x0_Callback(hObject, eventdata, handles)
plotgraph(handles);

function mu_Callback(hObject, eventdata, handles)
plotgraph(handles);

function position_Callback(hObject, eventdata, handles)
plotgraph(handles);

function velocity_Callback(hObject, eventdata, handles)
plotgraph(handles);

function minimumJerk_Callback(hObject, eventdata, handles)
plotgraph(handles);

function lognormal_Callback(hObject, eventdata, handles)
plotgraph(handles);

function sigma_Callback(hObject, eventdata, handles)
plotgraph(handles);

function sigma2_Callback(hObject, eventdata, handles)
plotgraph(handles);

function mu2_Callback(hObject, eventdata, handles)
plotgraph(handles);

function D2_Callback(hObject, eventdata, handles)
plotgraph(handles);

function A2_Callback(hObject, eventdata, handles)
plotgraph(handles);

function t02_Callback(hObject, eventdata, handles)
plotgraph(handles);

function x02_Callback(hObject, eventdata, handles)
plotgraph(handles);

function onesubmovement_Callback(hObject, eventdata, handles)
plotgraph(handles);

function twosubmovements_Callback(hObject, eventdata, handles)
plotgraph(handles);


function A_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function x0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function D_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function t0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sigma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function D2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function A2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function t02_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function x02_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function mu2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function sigma2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
