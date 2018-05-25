function [posEst,linVelEst,oriEst,driftEst,...
    posVar,linVelVar,oriVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,driftEst,...
%    posVar,linVelVar,oriVar,driftVar,estState] =
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConstants.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2018
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch
%
% Solutions coded by Guoxiang Zhou, Ruben Buitenhuis and Hong Fai Lau

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    % initial state mean
    posEst = [0,0]; % 1x2 matrix
    linVelEst = [0,0]; % 1x2 matrix
    oriEst = 0.5*(-estConst.RotationStartBound+estConst.RotationStartBound); % 1x1 matrix
    driftEst = 0.5*(-estConst.GyroDriftStartBound+estConst.GyroDriftStartBound); % 1x1 matrix
    
    % initial state variance
    posVar = [estConst.StartRadiusBound^2/4,estConst.StartRadiusBound^2/4]; % 1x2 matrix / calculated from joint pdf: p(x,y)=1/pi*R_0^2, and marginalizing to find p(x) and p(y)
    linVelVar = [0,0]; % 1x2 matrix / known
    oriVar = 1/12*(estConst.RotationStartBound-(-estConst.RotationStartBound))^2; % 1x1 matrix / uniform distribution
    driftVar = 1/12*(estConst.GyroDriftStartBound-(-estConst.GyroDriftStartBound))^2; % 1x1 matrix / uniform distribution
    
    % estimator variance init (initial posterior variance)
    estState.Pm = [posVar,linVelVar,oriVar,driftVar]'; % keep the diagonal elements, column vector
    % estimator state
    estState.xm = [posEst,linVelVar,oriVar,driftVar]';
    % time of last update
    estState.tm = tm;
    return;
end

%% Estimator iteration.

% get time since last estimator update
dt = tm-estState.tm;
estState.tm = tm; % update measurement time

%% prior update
Qc = diag([estConst.DragNoise,estConst.RudderNoise,estConst.GyroDriftNoise]); % process model noise variance matrix

% ODE solver for Ricatti Equation, using helper function
[xp, Pp] = priorUpdate(Qc, actuate, dt, estConst, estState);

%% measurement update
% H and M as defined in Lecture 8
H = [(xp(1)-estConst.pos_radioA(1))/sqrt((xp(1)-estConst.pos_radioA(1))^2+(xp(2)-estConst.pos_radioA(2))^2) (xp(2)-estConst.pos_radioA(2))/sqrt((xp(1)-estConst.pos_radioA(1))^2+(xp(2)-estConst.pos_radioA(2))^2) 0 0 0 0;
    (xp(1)-estConst.pos_radioB(1))/sqrt((xp(1)-estConst.pos_radioB(1))^2+(xp(2)-estConst.pos_radioB(2))^2) (xp(2)-estConst.pos_radioB(2))/sqrt((xp(1)-estConst.pos_radioB(1))^2+(xp(2)-estConst.pos_radioB(2))^2) 0 0 0 0;
    (xp(1)-estConst.pos_radioC(1))/sqrt((xp(1)-estConst.pos_radioC(1))^2+(xp(2)-estConst.pos_radioC(2))^2) (xp(2)-estConst.pos_radioC(2))/sqrt((xp(1)-estConst.pos_radioC(1))^2+(xp(2)-estConst.pos_radioC(2))^2) 0 0 0 0;
    0                                                                                                      0                                                                                                      0 0 1 1;
    0                                                                                                      0                                                                                                      0 0 1 0];
M = eye(5);

% all sensor noises are mutually independent, so R is diagonal:
R = diag([estConst.DistNoiseA,estConst.DistNoiseB,estConst.DistNoiseC,estConst.GyroNoise,estConst.CompassNoise]);

% measurement model
hk = [sqrt((xp(1)-estConst.pos_radioA(1))^2+(xp(2)-estConst.pos_radioA(2))^2);
    sqrt((xp(1)-estConst.pos_radioB(1))^2+(xp(2)-estConst.pos_radioB(2))^2);
    sqrt((xp(1)-estConst.pos_radioC(1))^2+(xp(2)-estConst.pos_radioC(2))^2);
    xp(5)+xp(6);
    xp(5)];

% check if measurement for radio C is missing, if so remove all
% corresponding entries in matrices.
if sum(isinf(sense))+sum(isnan(sense)) > 0
    H(3,:) = [];
    M = eye(4);
    R = diag([estConst.DistNoiseA,estConst.DistNoiseB,estConst.GyroNoise,estConst.CompassNoise]);
    hk(3) = [];
    sense(3) = [];
end

% HEKF measurement update equations, referencing Lecture 8
[estState.xm, estState.Pm] = measurementUpdate(xp,Pp,H,M,R,sense,hk)

%% Output quantities
posEst = estState.xm(1:2);
linVelEst = estState.xm(3:4);
oriEst = estState.xm(5);
driftEst = estState.xm(6);

posVar = estState.Pm(1:2);
linVelVar =  estState.Pm(3:4);
oriVar =  estState.Pm(5);
driftVar = estState.Pm(6);

end

%% functions
function dq = odefunSystem(q,Qc,actuate,sizeX,estConst)
% Helper function for solving the ODE, allows to solve for the system
% states and variance matrix P simultaneously. To allow ODE45 to work,
% matrix dP is reshaped as a vector and added below dx, which form the
% output dq. The input q is a vector with the states x and matrix P in
% vector form.

% reshape P as a matrix from q, and extract x from q.
% P = reshape(q(sizeX+1:end),sqrt(length(q(sizeX+1:end))),sqrt(length(q(sizeX+1:end))));
P = diag(q(sizeX+1:end));

x = q(1:sizeX);

% process model
dx = [ x(3);
    x(4);
    cos(x(5))*(tanh(actuate(1))-estConst.dragCoefficient*(x(3)^2+x(4)^2)*(1+0));
    sin(x(5))*(tanh(actuate(1))-estConst.dragCoefficient*(x(3)^2+x(4)^2)*(1+0));
    estConst.rudderCoefficient*actuate(2)*(1+0);
    0 ];

% define time-varying matrices A and L
A =  [ 0 0 1                                          0                                          0                                                                      0;
    0 0 0                                          1                                          0                                                                      0;
    0 0 -2*x(3)*estConst.dragCoefficient*cos(x(5)) -2*x(4)*estConst.dragCoefficient*cos(x(5)) -sin(x(5))*(tanh(actuate(1))-estConst.dragCoefficient*(x(3)^2+x(4)^2)) 0;
    0 0 -2*x(3)*estConst.dragCoefficient*sin(x(5)) -2*x(4)*estConst.dragCoefficient*sin(x(5)) cos(x(5))*(tanh(actuate(1))-estConst.dragCoefficient*(x(3)^2+x(4)^2))  0;
    0 0 0                                          0                                          0                                                                      0;
    0 0 0                                          0                                          0                                                                      0];
L =  [ 0                                                     0                                     0;
    0                                                     0                                     0;
    cos(x(5))*(-estConst.dragCoefficient*(x(3)^2+x(4)^2)) 0                                     0;
    sin(x(5))*(-estConst.dragCoefficient*(x(3)^2+x(4)^2)) 0                                     0;
    0                                                     estConst.rudderCoefficient*actuate(2) 0;
    0                                                     0                                     1];

% Ricatti Equation for P:
dP = A*P+P*A'+L*Qc*L';

% Build output vector dq
dq = [dx;diag(dP)];

end

function [xp, Pp]= priorUpdate(Qc, actuate, dt, estConst, estState)
% Qc is the process model noise variance matrix
% actuate is the control input
% dt: get time since last estimator update
% estState: previous posterior state
[~,qOut] = ode45(@(t,q)odefunSystem(q,Qc,actuate,6,estConst),[estState.tm-dt,estState.tm],[estState.xm;estState.Pm]);
qOut = qOut';
qOut = qOut(:,end);
xp = qOut(1:6); qOut(1:6) = [];
Pp = diag(qOut);
end

function [xm, Pm] = measurementUpdate(xp,Pp,H,M,R,sense,hk)
% xp: prior state
% xp: prior state variance matrix
% H,M,R,hk : as defined in lecture 8
% sence: measurements at current time
K = Pp*H'/(H*Pp*H' + M*R*M');
xm = xp + K*(sense'-hk);
Pm_maxtrix = (eye(size(Pp))-K*H)*Pp;
Pm = diag(Pm_maxtrix);
end
