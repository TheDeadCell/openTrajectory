%1-DOF Trajectory Simulation Model for OpenMotor Data
%Authored by Deouskar. A
%Salient Features: ISA Atmosphere with altitude varying g, speed of sound
%Salient Features: Mach-dependent Cd, thrust and propellant mass read from an OpenMotor CSV.
%Salient Features: Automatic mass handling for multi-grain motors.
clear;
clc;
close all;


%User Input Engine:
fprintf('=========== 1-DOF TRAJECTORY SIMULATOR ===========\n\n');
dryMass=askPositive('Dry Mass of the Rocket (i.e. mass of rocket w/o prop): ');
BodyDim=askPositive('Diameter of the body tube in mm: ');
crossArea = pi * (BodyDim/1000)^2 / 4;                  % [m^2] reference area
Cd0     = askPositive('Subsonic drag coefficient Cd (typical 0.3-0.5): ');
h0      = askNumber('Launch site elevation above sea level [m] (Enter = 0): ', 0);

%----------------------------------------------------------------------------------------------------------------------------------------

%File Opener: 
[fn, fp]=uigetfile({'*.csv','OpenMotor CSV export (*.csv)'}, 'Select the OpenMotor CSV file');
if isequal(fn, 0), error('No file selected.'); end
Tbl   = readtable(fullfile(fp, fn), 'VariableNamingRule', 'preserve');
names = Tbl.Properties.VariableNames;


%Load Operations
iT  = pickColumn(names, {'time'},           {}, 'TIME');
iF  = pickColumn(names, {'force','thrust'}, {}, 'THRUST');
iMs = pickMassColumns(names);       % ALL propellant-mass columns (one per grain)

fprintf('\nColumns detected:\n');
fprintf('  time   : %s\n', names{iT});
fprintf('  thrust : %s\n', names{iF});
for k = 1:numel(iMs)
    fprintf('  mass %d : %s\n', k, names{iMs(k)});
end
fprintf('  -> %d grain mass column(s) will be summed.\n', numel(iMs));

% Units are read from each header, e.g. "Thrust(N)" or "Propellant Mass(G1;kg)"
t_csv = double(Tbl{:, iT}) * unitScale(names{iT}, 'time');
F_csv = double(Tbl{:, iF}) * unitScale(names{iF}, 'force');
m_grain = zeros(height(Tbl), numel(iMs));
for k = 1:numel(iMs)
    m_grain(:, k) = double(Tbl{:, iMs(k)}) * unitScale(names{iMs(k)}, 'mass');
end
m_csv = sum(m_grain, 2);            % TOTAL propellant mass remaining (all grains)

% Clean: remove NaNs, sort, make time strictly increasing
good = ~(isnan(t_csv) | isnan(F_csv) | isnan(m_csv));
[t_csv, order] = sort(t_csv(good));
F_csv = F_csv(good); F_csv = F_csv(order);
m_csv = m_csv(good); m_csv = m_csv(order);
[t_csv, iu] = unique(t_csv);  F_csv = F_csv(iu);  m_csv = m_csv(iu);

burnoutTime = t_csv(end);                          % [s] burn ends with the last CSV row
TotalI      = trapz(t_csv, F_csv);                 % [N s] total impulse
F_avg       = TotalI / (t_csv(end) - t_csv(1));
massProp0   = m_csv(1);
iStart      = find(F_csv > 0, 1);                  % first non-zero thrust sample (OpenMotor logs 0 N at t = 0)
if isempty(iStart), error('Thrust column has no positive values.'); end
TWR0        = F_csv(iStart) / ((dryMass + massProp0) * 9.80665);

fprintf('\n--- Motor data summary ---\n');
fprintf('Burn time:        %.3f s  (%d data points, dt ~ %.3f s)\n', burnoutTime, numel(t_csv), median(diff(t_csv)));
fprintf('Peak thrust:      %.1f N\n', max(F_csv));
fprintf('Average thrust:   %.1f N\n', F_avg);
fprintf('Total impulse:    %.1f N s\n', TotalI);
fprintf('Propellant mass:  %.3f kg initial, %.3f kg at last row\n', massProp0, m_csv(end));
for k = 1:numel(iMs)
    fprintf('   grain %d: %.3f kg initial\n', k, m_grain(1, k));
end
fprintf('Liftoff mass:     %.3f kg   |  Initial thrust-to-weight: %.2f\n', dryMass + massProp0, TWR0);
if TWR0 < 1
    warning('Initial thrust-to-weight < 1: the rocket will not lift off.');
end

% Simulation Parameters Defination-----------------------------------------

dt = 0.002;                 % [s] integration time step
maxTime=300;                % [s] maximum simulated time
time = 0 : dt : maxTime;
N=numel(time);

velocity=zeros(1,N);
altitude=zeros(1,N);        % altitude above the launch site (AGL)
acceleration=zeros(1,N);
thrustProfile=zeros(1,N);
massProfile=zeros(1,N);
CdProfile=zeros(1,N);

apogeeAlt=0;
timeApogee=0;
apogeeReached=false;
launched=false;
maxMach=0;
maxAccel=0;

vBurnout=NaN;
hBurnout=NaN;
warnedAlt=false;
last=N;

%Numerical Integration Method --------------------------------------------

for i = 1:(N-1)
    t=time(i);
    vCurr=velocity(i);
    alt=altitude(i);

    % Apogee: the rocket has launched and its vertical speed is back to zero or below
    if launched && vCurr <=0
        apogeeAlt=alt;
        timeApogee=t;
        apogeeReached=true;
        last=i;
        break;
    end
    if vCurr>0, launched=true; end

    %the ISA Model with altitude dependent speed of sound and g----------------

    h_asl=h0+max(alt,0);
    if h_asl>32000 && ~warnedAlt
        warning('Altitude exceeds 32Km, the programmed ISA model may not be valid anymore');
        warnedAlt=true;
    end
    [rho, aSound, g]=atmosphereModel(h_asl); 


    %Drag Coefficient 

    Mach=abs(vCurr)/aSound; 
    maxMach=max(maxMach, Mach);
    Cd= dragCoefficient(Mach, Cd0);

    %Thrust and propellent mass from openMotor CSV

    if t<burnoutTime
        T=interp1(t_csv, F_csv, t, 'linear', 0);
        T=max(0,T);
    else
        T=0;
    end

    tq     = min(max(t, t_csv(1)), t_csv(end));          % clamp: hold last value after burnout
    m_prop = max(0, interp1(t_csv, m_csv, tq, 'linear'));
    m_total = dryMass + m_prop;

    % --- 1-D Forces --- %
    dragF=0.5 * rho * vCurr * abs(vCurr) * crossArea * Cd;
    netF  = T - dragF - m_total * g;
    a_curr = netF/ m_total;

    if alt <= 0 && vCurr <= 0 && netF < 0
        a_curr = 0;                                      % rocket rests on the launchpad
    end

    thrustProfile(i) = T;   massProfile(i) = m_total;
    CdProfile(i)     = Cd;  acceleration(i) = a_curr;
    maxAccel = max(maxAccel, (T - dragF)/m_total);       % g-load felt by the vehicle = (thrust - drag)/mass


    % --- Integration (trapezoidal position update) ---
    velocity(i+1) = vCurr + a_curr * dt;
    altitude(i+1) = alt + 0.5 * (vCurr + velocity(i+1)) * dt;

    if isnan(vBurnout) && t >= burnoutTime
        vBurnout = vCurr;  hBurnout = alt;
    end
end

if ~launched
    error('Rocket has not left the ground. Please check masses and compute motor thrust correctly.');
elseif ~apogeeReached
    warning('Not reached Apogee within %.0f s. Increase maximum simulation time (maxTime).', maxTime);
    apogeeAlt = max(altitude);  timeApogee = time(last);
end

time = time(1:last);            velocity = velocity(1:last);
altitude = altitude(1:last);    acceleration = acceleration(1:last);
thrustProfile = thrustProfile(1:last);
massProfile   = massProfile(1:last);
CdProfile     = CdProfile(1:last);


%% 5. RESULTS
fprintf('\n================ TRAJECTORY RESULTS ================\n');
fprintf('Duration of Motor Burn:    %.2f s\n', burnoutTime);
fprintf('Burnout velocity:       %.1f m/s at %.0f m AGL\n', vBurnout, hBurnout);
fprintf('Max velocity:           %.1f m/s (Mach %.2f)\n', max(velocity), maxMach);
fprintf('Max acceleration:       %.1f m/s^2 (%.1f g)\n', maxAccel, maxAccel/9.80665);
fprintf('Time to apogee:         %.2f s\n', timeApogee);
fprintf('PREDICTED APOGEE (AGL): %.1f m (%.2f km)\n', apogeeAlt, apogeeAlt/1000);
fprintf('====================================================\n');

%% Caution!: The below mentioned code related to plotting the graph has been generated using generative Artifical Intelligence!
figure('Color', 'w', 'Position', [100 100 800 800]);

subplot(3,1,1);
yyaxis left
plot(time, altitude/1000, 'b-', 'LineWidth', 2.5); hold on;
plot(timeApogee, apogeeAlt/1000, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(timeApogee, apogeeAlt/1000, sprintf('  Apogee: %.2f km', apogeeAlt/1000), ...
    'Color', 'r', 'FontWeight', 'bold', 'VerticalAlignment', 'top');
ylabel('Altitude AGL (km)');
yyaxis right
plot(time, velocity, 'r--', 'LineWidth', 2);
ylabel('Velocity (m/s)');
grid on; title('Flight trajectory to apogee');

subplot(3,1,2);
yyaxis left
plot(time, thrustProfile, 'k-', 'LineWidth', 2);
ylabel('Thrust (N)');
yyaxis right
plot(time, massProfile, 'm-', 'LineWidth', 2);
ylabel('Total mass (kg)');
grid on; title('Motor thrust (OpenMotor CSV) and vehicle mass');
xlim([0, min(time(end), burnoutTime*3)]);

subplot(3,1,3);
plot(time, CdProfile, 'k-', 'LineWidth', 2);
grid on; ylabel('Cd'); xlabel('Time (s)');
title('Drag coefficient');


% Mach number along the trajectory (needed for the altitude-velocity graph)
machNo = zeros(size(altitude));
for i = 1:numel(altitude)
    [~, aLocal, ~] = atmosphereModel(h0 + max(altitude(i),0));
    machNo(i) = velocity(i)/aLocal;
end

% --- Graph 1: Altitude vs Time ---
figure('Color', 'w');
plot(time, altitude/1000, 'b-', 'LineWidth', 2.5); hold on;
plot(timeApogee, apogeeAlt/1000, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(timeApogee, apogeeAlt/1000, sprintf('Apogee: %.2f km  ', apogeeAlt/1000), ...
    'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
xline(burnoutTime, 'k--', 'Burnout');
grid on; xlabel('Time (s)'); ylabel('Altitude AGL (km)');
title('Altitude vs Time');

% --- Graph 2: Velocity vs Time ---
figure('Color', 'w');
plot(time, velocity, 'r-', 'LineWidth', 2.5); hold on;
[vMax, iMax] = max(velocity);
plot(time(iMax), vMax, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
text(time(iMax), vMax, sprintf('  Max: %.0f m/s (Mach %.2f)', vMax, machNo(iMax)), ...
    'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
xline(burnoutTime, 'k--', 'Burnout');
grid on; xlabel('Time (s)'); ylabel('Velocity (m/s)');
title('Velocity vs Time');

% --- Graph 3: Altitude vs Velocity, with Mach number markers ---
% (Mach is not proportional to velocity because the speed of sound changes with
%  altitude, so Mach is marked on the curve instead of on a second axis.)
figure('Color', 'w');
plot(velocity, altitude/1000, 'b-', 'LineWidth', 2.5); hold on;
for M = 0.5 : 0.5 : maxMach
    iUp   = find(machNo >= M, 1);              % first time the rocket reaches this Mach (accelerating)
    iDown = find(machNo >= M, 1, 'last');      % last time it is above this Mach (decelerating)
    for iM = unique([iUp, iDown])
        plot(velocity(iM), altitude(iM)/1000, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
        text(velocity(iM), altitude(iM)/1000, sprintf('  Mach %.1f', M), 'FontSize', 8);
    end
end
grid on; xlabel('Velocity (m/s)'); ylabel('Altitude AGL (km)');
title('Altitude vs Velocity (red markers = Mach number)');


%%The AI generated code ends here 

%%Physics Governing Functions (do not change the position of this block of code in the file! )

function [rho, a, g]=atmosphereModel(h)
% ISA atmosphere, valid 0-32 km (geopotential altitude), plus g(h).

R   = 287.053;   gam = 1.4;   g0 = 9.80665;   Re = 6371000;
H = Re*h/(Re + h);

if H < 11000
    L = -0.0065;
    T = 288.15 + L*H;
    p = 101325*(T/288.15)^(-g0/(R*L));
elseif H < 20000
    T = 216.65;
    p = 22632.06*exp(-g0*(H - 11000)/(R*T));
else
    L = 0.001;
    T = 216.65 + L*(H - 20000);
    p = 5474.889*(T/216.65)^(-g0/(R*L));
end
rho = p/(R*T);
a   = sqrt(gam*R*T);
g   = g0*(Re/(Re + h))^2;
end

function Cd= dragCoefficient(M, Cd0)
%Transonic Drag Coefficient Model-automatically approximates and scales the
%user defined Cd to match appropriate regime specific trends
s       = Cd0 / 0.35;
    Cd_peak = 0.75;
    Cd_13   = 0.55 / sqrt(1.3^2 - 1);
    if M < 0.8
        base = 0.35;
    elseif M < 1.05
        base = 0.35 + 0.40*((M - 0.8)/0.25)^2;
    elseif M <= 1.3
        base = Cd_peak + (Cd_13 - Cd_peak)*((M - 1.05)/0.25);
    else
        base = max(0.40, 0.55/sqrt(M^2 - 1));
    end
    Cd = s * base;
end

function idx=pickColumn(names, mustHave, mustNot, label)
%% Caution!: The below mentioned code related to reading the CSV files has been generated using generative Artifical Intelligence!

% Auto-detect a CSV column by header keywords; ask the user if not found.
ln = lower(names);
ok = false(size(ln));
for k = 1:numel(ln)
    has = any(cellfun(@(s) contains(ln{k}, s), mustHave));
    bad = ~isempty(mustNot) && any(cellfun(@(s) contains(ln{k}, s), mustNot));
    ok(k) = has && ~bad;
end
idx = find(ok, 1);
if isempty(idx)
    fprintf('\nCould not auto-detect the %s column. Columns found:\n', label);
    for k = 1:numel(names), fprintf('  %d: %s\n', k, names{k}); end
    idx = round(askNumber(sprintf('Enter column number for %s: ', label), NaN));
end
end

function x = askNumber(prompt, default)
while true
    s = strtrim(input(prompt, 's'));
    if isempty(s) && ~isnan(default), x = default; return; end
    x = str2double(s);
    if ~isnan(x), return; end
    fprintf('  Please enter a number.\n');
end
end

function x = askPositive(prompt)
while true
    x = askNumber(prompt, NaN);
    if x > 0, return; end
    fprintf('  Value must be greater than 0.\n');
end
end

function idx = pickMassColumns(names)
% Find every propellant-mass column (one per grain in OpenMotor exports).
% If a "total" mass column exists it is used alone, to avoid double counting.
ln     = lower(names);
isMass = contains(ln, 'mass') & ~contains(ln, 'flow') & ~contains(ln, 'flux');
isTot  = isMass & contains(ln, 'total');
if any(isTot)
    idx = find(isTot, 1);
else
    idx = find(isMass);
end
if isempty(idx)
    fprintf('\nCould not auto-detect mass columns. Columns found:\n');
    for k = 1:numel(names), fprintf('  %d: %s\n', k, names{k}); end
    s = input('Enter mass column number(s), e.g. [3 4]: ');
    idx = round(s(:)');
end
end

function sc = unitScale(header, kind)
% Read the unit from the header text in the last parentheses
% ("Thrust(N)", "Propellant Mass(G1;kg)") and return the factor to SI.
tok = regexp(header, '\(([^)]*)\)\s*$', 'tokens', 'once');
u = '';
if ~isempty(tok)
    parts = strsplit(tok{1}, ';');
    u = lower(strtrim(parts{end}));
end
switch kind
    case 'time',  tbl = {'s',1; 'ms',1e-3; 'min',60};
    case 'force', tbl = {'n',1; 'kn',1e3; 'lbf',4.4482216};
    case 'mass',  tbl = {'kg',1; 'g',1e-3; 'lb',0.45359237; 'lbs',0.45359237; 'oz',0.028349523};
end
k = find(strcmp(u, tbl(:,1)), 1);
if isempty(k)
    warning('Unit "%s" in header "%s" not recognised; assuming SI.', u, header);
    sc = 1;
else
    sc = tbl{k, 2};
end
end
