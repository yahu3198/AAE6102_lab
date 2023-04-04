%##########################################################################
%                           AAE6102 Lab
%                        Yang Hu 22073277R
%##########################################################################

%% This lab focus on calculating user position and clock bias based on ephemeris data and pseudorange
%% input data
clear;clc;
% Load data from files
% addpath('/home/fyp/GNSS_lab/Lab/Data');
load rcvr.dat
load eph.dat

eph = sortrows(eph,2);
rcvr = sortrows(rcvr,2);

% Input parameters
current_pos = [-2694685.473;...             % Initial user position with clock bias
               -4293642.366;...
               3857878.924;...
               0];
thre_dx = 1e-4;                             % Threshold to stop the iteration

ref_pos = [-2700400; -4292560; 3855270];    % Reference user position

c = 299792458;                              % Speed of light (m/s)
wedot = 7.2921151467e-5;                    % WGS 84 value of earth's rotation rate (r/s)
mu = 3.986005e+14;                          % WGS 84 value of earth's universal gravitation constant (m^3/s^2)
F = -4.442807633e-10;                       % Relativistic correction term constant


%% calculate satellites position and clock bias (ICD file Table 20-IV)
for i = 1:size(eph,1)
    rcvr_tow    = eph(i,1);     % Receiver time of week (s)
    svid        = eph(i,2);     % Satellite PRN number (1-32)
    toc         = eph(i,3);     % Reference time of clock parameters (s)
    toe         = eph(i,4);     % Reference time of ephemeris parameters (s)
    af0         = eph(i,5);     % Clock correction coefficient - group delay (s)
    af1         = eph(i,6);     % Clock correction coefficient (s/s)
    af2         = eph(i,7);     % Clock correction coefficient (s/s/s)
    ura         = eph(i,8);     % User range accuracy (m)
    e           = eph(i,9);     % Eccentricity (-)
    sqrta       = eph(i,10);	% Square root of semi-major axis a (m^1/2)
    dn          = eph(i,11);	% Mean motion correction (r/s)
    m0          = eph(i,12);	% Mean anomaly at reference time (r)
    w           = eph(i,13);	% Argument of perigee (r)
    omg0        = eph(i,14);	% Right ascension (r)
    i0          = eph(i,15);    % Inclination angle at reference time (r)
    odot        = eph(i,16);    % Rate of right ascension (r/s)
    idot        = eph(i,17);    % Rate of inclination angle (r/s)
    cus         = eph(i,18);    % Argument of latitude correction, sine (r)
    cuc         = eph(i,19);    % Argument of latitude correction, cosine (r)
    cis         = eph(i,20);    % Inclination correction, sine (r)
    cic         = eph(i,21);    % Inclination correction, cosine (r)
    crs         = eph(i,22);    % Radius correction, sine (r)
    crc         = eph(i,23);    % Radius correction, cosine (r)
    iod         = eph(i,24);    % Issue of data number
    
    a = sqrta^2;                % Semi-major axis
    n0 = sqrt(mu/a^3);          % Computed mean motion (r/s)
    t = rcvr_tow - rcvr(i,3)/c; % Satellite signal transmition time
    tk = t - toe;               % Time from ephemeris reference epoch
    
    % Account for beginning or end of week crossovers
    if tk > 302400
        tk = tk - 604800;
    elseif tk < -30240
        tk = tk + 604800;
    end
    
    n = n0+dn;                  % Corrected mean motion
    mk = m0+n*tk;               % Mean anomaly
    
    % Kepler's equation for Eccentric Anomaly (solved by iteration)
    iter = 0;
    ek = mk;
    ek0 = ek;
    while iter<12 || abs(ek-ek0)>1e-12
        ek0 = ek;
        ek = ek - (ek - e*sin(ek) - mk)/(1-e*cos(ek));
        iter = iter + 1;
    end
    
    vk = atan2((sqrt(1-e^2)*sin(ek)/(1-e*cos(ek))),((cos(ek)-e)/(1-e*cos(ek))));    % True anomaly
    
    phik = vk + w;                              % Argument of latitude
    
    % Second harmonic perturbations
    duk = cus*sin(2*phik) + cuc*cos(2*phik);    % Argument of Latitude Correction
    drk = crs*sin(2*phik) + crc*cos(2*phik);    % Radius Correction
    dik = cis*sin(2*phik) + cic*cos(2*phik);    % Inclination Correction
    
    uk = phik + duk;                            % Corrected argument of latitude
    rk = a*(1 - e*cos(ek)) + drk;               % Corrected Radius
    ik = i0 + dik + idot*tk;                    % Corrected Inclination
    xkp = rk*cos(uk);                           % X position in orbital plane
    ykp = rk*sin(uk);                           % Y position in orbital plane
    omgk = omg0 + (odot - wedot)*tk -wedot*toe; % Corrected longitude of ascending node
    
    % Earth rotation correction compensation
    xk = xkp*cos(omgk) - ykp*cos(ik)*sin(omgk); % Earth-fixed X coordinates
    yk = xkp*sin(omgk) + ykp*cos(ik)*cos(omgk); % Earth-fixed Y coordinates
    zk = ykp*sin(ik);                           % Earth-fixed Z coordinates
    
    % Calculate satellite clock bias (s)
    dtsv = af0 + af1*(t - toc) + af2*(t - toc)^2 + F*e*sqrta*sin(ek);
    
    sat_data(i,:) = [xk yk zk dtsv];
end


%% get user position and clock bias, and print out
% Correct the pseudorange
sat_data(:,4) = rcvr(:,3) + c .* sat_data(:,4);

% Use least square to estimate user position and clock bias
dx = [-5710, 1080, -2610, 519450];  % Initial iteration
iter = 0;                           % Number of iterations

disp('#######################################################################################')
disp('--------------------------------Iteration starts---------------------------------------')
fprintf('Initial user position in earth-fixed coordinates and user clock bias\n');
fprintf('X coordinate: %f (m)\n',current_pos(1));
fprintf('Y coordinate: %f (m)\n',current_pos(2));
fprintf('Z coordinate: %f (m)\n',current_pos(3));
fprintf('User clock bias: %f (s)\n',current_pos(4)/c);
fprintf('\n');

while norm(dx(1:3)) > thre_dx
    dx = estimate_user_pos(sat_data,current_pos);
    current_pos = current_pos + dx;   % Update approximation
    iter = iter+1;
    disp('---------------------------------------------------------------------------------------')
    fprintf('Iteration No.%d result:\n',iter);
    fprintf('Updated user position in earth-fixed coordinates and user clock bias\n');
    fprintf('X coordinate: %f (m)\n',current_pos(1));
    fprintf('Y coordinate: %f (m)\n',current_pos(2));
    fprintf('Z coordinate: %f (m)\n',current_pos(3));
    fprintf('User clock bias: %f (s) / %f (m)\n',current_pos(4)/c, current_pos(4));
    fprintf('\n');
end

% Print the results
disp('-------------------------Threshold reached: iteration ends-----------------------------')
disp('#######################################################################################')



%% function of estimating user position and clock bias
function [dx] = estimate_user_pos(sat_data,current_pos)
    
    % Input parameters
    sat_x = sat_data(:,1);      % Satellite X position
    sat_y = sat_data(:,2);      % Satellite Y position
    sat_z = sat_data(:,3);      % Satellite Z position
    pr = sat_data(:,4);         % Measured pseudorange
    x0 = current_pos(1);        % Initial X user position
    y0 = current_pos(2);        % Initial Y user position
    z0 = current_pos(3);        % Initial Z user position
    b0 = current_pos(4);        % Initial user clock bias (m)
    
    rho = pr;                   % Current (measuren) pseudorange
    r = sqrt((sat_x - x0).^2 + (sat_y - y0).^2 + (sat_z - z0).^2); % Geometric range between user and satellite
                           
    rho_hat = r - b0;           % Last (approximated) pseudorange
    delta_rho = rho_hat - rho;
    
    % Calculate H matrix
    ax = (sat_x - x0)./r;
    ay = (sat_y - y0)./r;
    az = (sat_z - z0)./r;
    H = [ax ay az ones(size(sat_data,1),1)];
    
    % Calculate dx
    dx = inv(H'*H)*H'*delta_rho;

end