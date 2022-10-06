% Bloch simulation of the slice-profile or time-course of the magnetisation
% resulting from a hamming-windowed sinc RF pulse and slice select gradient 

%% User inputs

% Chose what to simulate with sim_type
% 1. slice: Simulate the slice profile immediately after the slice selection RF pulse and gradient
% 2. time: Simulate mx, my and mz throughout the slice-selection process
% 3. freq_pos: Simulate the slice profile over a range of off-resonance frequencies
sim_type = 'freq_pos';  %slice, time, freq_pos

gamma = 2.6753E8;       % Proton gyromagnetic ratio (rad s-1 T-1)

% Properties of the slice selection pulse
flip_rf = pi / 2;       % Flip-angle (radians)
pw_rf = 3.2E-3;         % Duration (s)
tbw = 8;                % Time-bandwidth product 
bw_rf = tbw / pw_rf;    % Bandwidth (Hz)
alpha_rf = 0.46;        % Sinc apodisation 0=none, 0.46=hamming, 0.5=hanning
n_rf = 800;             % Number of points
dt = pw_rf / n_rf;      % Time interval (s)
slthick = 5E-3;         % Slice thickness (m)

t1=1;                   % Longitudinal relaxation (s)
t2=0.2;                 % Transverse relaxation (s)


%% Generate RF pulse and gradient waveform
b1 = msinc(n_rf, tbw / 4, alpha_rf);             % Pulse waveform (unitless)
a_rf = 10000 * flip_rf / (sum(b1) * dt * gamma); % RF amplitude (G)
b1 = [a_rf .* b1 zeros(1, n_rf / 2)];            % RF pulse waveform (G)

a_gz = 100 * 2 * pi * bw_rf / (gamma * slthick); % (G/cm)
gz=[ones(1, n_rf) -ones(1, n_rf / 2)] .* a_gz;   % Gradient waveform (G/cm)
g = zeros(n_rf + n_rf / 2, 3);
g(:, 3) = gz';

t = [1 : length(b1)] * dt * 1000; % time (ms)

figure(1);
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
if strcmp(sim_type, 'time')
    r = zeros(1, 3); % Position vector
    f = 0; % Vector of offset frequencies
    mode = 2; % Simulation mode: From start recording all time points
elseif strcmp(sim_type, 'slice') 
    n_z = 1000; %The number of z-positions to simulate
    z = linspace(-200 .* slthick, 200 .* slthick, n_z);
    r = [zeros(n_z, 1),zeros(n_z, 1), z'];
    f = 0;
    mode = 0; % Simulation mode: From start, but only record the end timepoint
elseif strcmp(sim_type, 'freq_pos')
    n_z = 100; %The number of z-positions to simulate
    z = linspace(-200 .* slthick, 200 .* slthick, n_z);
    r = [zeros(n_z, 1), zeros(n_z, 1), z'];
    n_freq = 100;
    f = linspace(-250, 250, n_freq);
    mode = 0; %Simulation mode: From start, but only record the end timepoint
end


%Run bloch simulation
[mx, my, mz] = bloch(b1, g, dt, t1, t2, f, r, mode);
mxy = mx + 1i * my;

figure(2);
if strcmp(sim_type, 'time')
    
    yyaxis left
    plot(t, my);
    xlabel('Time (ms)');
    ylabel('M_y');
        
    yyaxis right;
    plot(t, mz);
    xlabel('Time (ms)');
    ylabel('M_z');
    
elseif strcmp(sim_type, 'slice')
    
    yyaxis left
    plot(z, abs(mxy));
    xlabel('Position (cm)');
    ylabel('|M_{xy}|');
    
    yyaxis right
    plot(z, angle(mxy));
    xlabel('Position (cm)');
    ylabel('Phase M_{xy} (rad)');
    
elseif strcmp(sim_type,'freq_pos')
    
    imagesc(f, z, abs(mxy(:,  :)));
    colormap('gray');
    colorbar;
    xlabel('Off Resonance Frequency (Hz)') 
    ylabel('Position (cm)');
    title('|M_{xy}|');
end

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


