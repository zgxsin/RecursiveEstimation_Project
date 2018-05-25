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

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    
    % initial state mean
    posEst = [0,0]; % 1x2 matrix
    linVelEst = [0,0]; % 1x2 matrix
    oriEst = 0.5*(-RotationStartBound+RotationStartBound); % 1x1 matrix
    driftEst = 0.5*(-GyroDriftStartBound+GyroDriftStartBound); % 1x1 matrix
    
    % initial state variance
    posVar = [StartRadiusBound^2/4,StartRadiusBound^2/4]; % 1x2 matrix
    linVelVar = [0,0]; % 1x2 matrix
    oriVar = 1/12*(RotationStartBound-(-RotationStartBound))^2; % 1x1 matrix
    driftVar = 1/12*(GyroDriftStartBound-(-GyroDriftStartBound))^2; % 1x1 matrix
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar,linVelVar,oriVar,driftVar]);
    % estimator state
    estState.xm = posEst;
    % time of last update
    estState.tm = tm;
    return;
end

%% Estimator iteration.

% get time since last estimator update
dt = tm-estState.tm;
estState.tm = tm; % update measurement update time

% prior update
% TODO: Define q, A(t), L(t).
%xp = q(xest,0,t) from (k-1)T<=t<=kT and x((k-1)T) estState.xm
%Pp = P(kT); Pdot(t) = A(t)*P(t)+P(t)*A'(t)+L(t)*Qc(t)*L'(t), for (k-1)T<=t<=kT and P((k-1)T) = estState.Pm;

% measurement update
K = Pp*H'/(H*Pp*H' + M*R*M');
estState.xm = xp + K*(sense-hk(xp,0)); % TODO: Define hk
estState.Pm = (eye(size(Pm))-K*H)*Pp;

% Get resulting estimates and variances
% Output quantities
posEst = estState.xm(1:2);
linVelEst = estState.xm(3:4);
oriEst = estState.xm(5);
driftEst = estState.xm(6);

posVar = [estState.Pm(1,1),estState.Pm(2,2)];
linVelVar = [estState.Pm(3,3),estState.Pm(4,4)];
oriVar = estState.Pm(5,5);
driftVar = estState.Pm(6,6);

end