% Bloch simulation of the slice-profile or time-course of the magnetisation
% resulting from a spectral-spatial RF pulse and gradient 

% See section 5.4.1 of Handbook of MRI Pulse Sequences for the underlying
% theory

%% User inputs

gamma = 2.6753E8;                               % Proton gyromagnetic ratio (rad s-1 T-1)
fcs = 250;                                      % Fat water chemical shift at 1.5 T (Hz) 
t1=1;                                           % Longitudinal relaxation (s)
t2=0.2;                                         % Transverse relaxation (s)

% Properties of the spatial RF sub-pulses
flip_rf = pi / 2;                               % Flip-angle (radians)
alpha_rf_spatial = 0.46;                        % Sinc apodisation (0.46=hamming)
n_rf_spatial = 800;                             % Number of points
slthick = 5E-3;                                 % Slice thickness (m)
tbw_spatial = 4;                                % Spatial time-bandwidth product
pw_rf_spatial = 1 / (4 * fcs);                  % Duration (s) see eq. 5.31
bw_rf_spatial = tbw_spatial / pw_rf_spatial;    % Spatial bandwidth (Hz)
dt = pw_rf_spatial / n_rf_spatial;              % Time interval (s)

% Properties of the spectral envelope
alpha_rf_spectral = 0.46;                       % Sinc apodisation
tbw_spectral = 4;                               % Spectral time-bandwidth product
bw_rf_spectral = fcs;                           % Spectral bandwidth (Hz)
pw_rf_spectral = tbw_spectral / bw_rf_spectral; % Duration (s)
n_rf_spectral = pw_rf_spectral / dt;            % Number of points

%% Generate RF pulse and gradient waveform

% Spatial sub-pulses 
b1_spatial_subpulse = msinc(n_rf_spatial, tbw_spatial / 4, alpha_rf_spatial);

% Replicate the correct number of spatial subpulses to fit in the spectral
% envelope
num_rf_pulses = ceil(pw_rf_spectral / pw_rf_spatial);
b1_spatial = repmat(b1_spatial_subpulse, 1, num_rf_pulses);

% Spectral envelope 
b1_spectral = msinc(n_rf_spectral, tbw_spectral / 4, alpha_rf_spectral);

% Gradient waveform
Gz = 2 * pi * bw_rf_spatial / (gamma .* slthick);  % Gradient amplitude (T/m)
gz = Gz .* ones(1, n_rf_spatial);
gz = [-gz gz];
gz = repmat(gz, 1, num_rf_pulses / 2);

g = zeros(length(gz), 3);
g(:, 3) = 100 .* gz'; % gradient waveform (G/cm)

% Normalise b1 using eqs. 5.33 and 5.34 to give desired flip angle
b1 = b1_spatial .* b1_spectral .* abs(gz);  % eq. 5.33
b1 = 10000 .* b1 .* flip_rf ./ (sum(b1) * dt * gamma); % eq. 5.34 (G)

figure(1);
t = [1 : length(b1)] * dt * 1000; % time (ms)
yyaxis left
plot(t, real(b1));
xlabel('Time (ms)');
ylabel('RF (G)');

yyaxis right
plot(t, g(:, 3));
xlabel('Time (ms)');
ylabel('Gradient (G/cm)');

%% Perform Bloch simulation

% Setup the parameters of the simulation
n_z = 200; % The number of z-positions to simulate
n_freq = 200;
z = linspace(-200 .* slthick, 200 .* slthick, n_z);
r = [zeros(n_z, 1), zeros(n_z, 1), z'];
f = linspace(-4 * fcs, 4 * fcs, n_freq);
mode = 0; % Simulation mode: From start, but only record the end timepoint

% Run bloch simulation
[mx, my, mz] = bloch(b1, g, dt, t1, t2, f, r, mode);
mxy = mx + 1i * my;

figure(2);
imagesc(f, z, abs(mxy(:,  :)));
colormap('gray');
colorbar;
xlabel('Off Resonance Frequency (Hz)') 
ylabel('Position (cm)');
title('|M_{xy}|');

%% Functions

function ms = msinc(n, m, alpha)
%  Return a windowed sinc of length n, with m sinc-cycles,
%  i.e. a time-bandwidth of 4*m, alpha 0=none, 0.46=hamming, 
%  0.5=hanning sinc apodisation
%
%  Originally written by John Pauly, 1992
%  (c) Board of Trustees, Leland Stanford Junior University
%
%  Modified to to include alpha by Stephen Wastling

x = [-n / 2 : (n - 1) / 2] / (n / 2);
snc = sin(m * 2 * pi * x + 0.00001) ./ (m * 2 * pi * x + 0.00001);
ms = snc .* ((1 - alpha) + alpha * cos(pi * x));
ms = ms * 4 * m / n;
end
