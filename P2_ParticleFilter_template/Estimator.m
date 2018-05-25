function [postParticles] = Estimator(prevPostParticles, sens, act, init)
% [postParticles] = Estimator(prevPostParticles, sens, act, init)
%
% The estimator function. The function will be called in two different
% modes: If init==1, the estimator is initialized. If init == 0, the
% estimator does an iteration for a single sample time interval Ts (KC.ts)
% using the previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurements and control inputs.
%
% You must edit this function.
%
% Inputs:
%   prevPostParticles   previous posterior particles at discrete time k-1,
%                       which corresponds to continuous time t = (k-1)*Ts
%                       The variable is a struct whose fields are arrays
%                       that correspond to the posterior particle states.
%                       The fields are: (N is number of particles)
%                       .x = (2xN) array with the x-locations (metres)
%                       .y = (2xN) array with the y-locations (metres)
%                       .h = (2xN) array with the headings (radians)
%                       The first row in the arrays corresponds to robot A.
%                       The second row corresponds to robot B.
%
%   sens                Sensor measurements at discrete time k (t = k*Ts),
%                       [4x1]-array, an Inf entry indicates no measurement
%                       of the corresponding sensor.
%                       sens(1): distance reported by sensor 1 (metres)
%                       sens(2): distance reported by sensor 2 (metres)
%                       sens(3): distance reported by sensor 3 (metres)
%                       sens(4): distance reported by sensor 4 (metres)
%
%   act                 Control inputs u at discrete time k-1, which are
%                       constant during a time interval Ts:
%                       u(t) = u(k-1) for (k-1)*Ts <= t < k*Ts
%                       [2x1]-array:
%                       act(1): velocity of robot A, u_A(k-1) (metres/second)
%                       act(2): velocity of robot B, u_B(k-1) (metres/second)
%
%   init                Boolean variable indicating wheter the estimator
%                       should be initialized (init = 1) or if a regular
%                       estimator update should be performed (init = 0).
%                       OPTIONAL ARGUMENT. By default, init = 0.
%
% Outputs:
%   postParticles       Posterior particles at discrete time k, which
%                       corresponds to the continuous time t = k*Ts.
%                       The variable is a struct whose fields are arrays
%                       that correspond to the posterior particle states.
%                       The fields are: (N is number of particles)
%                       .x = (2xN) array with the x-locations (metres)
%                       .y = (2xN) array with the y-locations (metres)
%                       .h = (2xN) array with the headings (radians)
%                       The first row in the arrays corresponds to robot A.
%                       The second row corresponds to robot B.
%
% Class:
% Recursive Estimation
% Spring 2018
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control

%% Solutions coded by Guoxiang Zhou, Hong Fai Lau and Ruben Buitenhuis

% Check if init argument was passed to estimator:
if(nargin < 4)
    % if not, set to default value:
    init = 0;
end

%% Mode 1: Initialization
% Set number of particles:
% obviously, you will need more particles than 10.
N = 5000; 
if floor(N) ~= N; error('N not integer'); end

if (init)
    % Do the initialization of your estimator here!
    % These particles are the posterior particles at discrete time k = 0
    % which will be fed into your estimator again at k = 1
    %Initializing Robot A at (2L,.) and Robot B at (0,.)
    postParticles.x = ones(2,N).*[2*KC.L;0]; 
    % y should be either 0 or KC.L
    postParticles.y = round(rand(2,N)).*KC.L;
    % Todo: I do not undersand
    postParticles.h = wrapToPi(rand(2,N).*pi/2 + ones(2,N).*[pi/2;0] + postParticles.y.*[-3*pi/2;-pi/2]./KC.L);
    return;
end % end init

%% Mode 2: Estimator iteration.
% If init = 0, we perform a regular update of the estimator.
% Implement your estimator here!
% Prior update/Prediction Step
for iloop = 1:N
    vsA = getVsRealization(KC.vsbar);
    vsB = getVsRealization(KC.vsbar);
    uA = act(1)*(1+vsA);
    uB = act(2)*(1+vsB);    
    
    priorParticles.x(1,iloop) = KC.ts*uA*cos(prevPostParticles.h(1,iloop))+prevPostParticles.x(1,iloop);
    priorParticles.x(2,iloop) = KC.ts*uB*cos(prevPostParticles.h(2,iloop))+prevPostParticles.x(2,iloop);
    % priorParticles.x(:,iloop) = KC.ts*[uA;uB].*cos(prevPostParticles.h(:,iloop))+prevPostParticles.x(:,iloop);
    
    priorParticles.y(1,iloop) = KC.ts*uA*sin(prevPostParticles.h(1,iloop))+prevPostParticles.y(1,iloop);
    priorParticles.y(2,iloop) = KC.ts*uB*sin(prevPostParticles.h(2,iloop))+prevPostParticles.y(2,iloop);
    
    % priorParticles.y(:,iloop) = KC.ts*[uA;uB].*cos(prevPostParticles.h(:,iloop))+prevPostParticles.y(:,iloop);
    
    priorParticles.h(:,iloop) = prevPostParticles.h(:,iloop);
    
    % Check if a bounce should happen
    for jloop = 1:2
        if priorParticles.x(jloop,iloop) < 0, [priorParticles.x(jloop,iloop),priorParticles.y(jloop,iloop),priorParticles.h(jloop,iloop)] = calculateWallBounce(priorParticles.x(jloop,iloop),priorParticles.y(jloop,iloop),priorParticles.h(jloop,iloop),-pi/2,0,[],KC.vbar); end
        if priorParticles.x(jloop,iloop) > KC.L*2, [priorParticles.x(jloop,iloop),priorParticles.y(jloop,iloop),priorParticles.h(jloop,iloop)] = calculateWallBounce(priorParticles.x(jloop,iloop),priorParticles.y(jloop,iloop),priorParticles.h(jloop,iloop),pi/2,KC.L*2,[],KC.vbar); end
        if priorParticles.y(jloop,iloop) < 0, [priorParticles.x(jloop,iloop),priorParticles.y(jloop,iloop),priorParticles.h(jloop,iloop)] = calculateWallBounce(priorParticles.x(jloop,iloop),priorParticles.y(jloop,iloop),priorParticles.h(jloop,iloop),pi,[],0,KC.vbar); end
        if priorParticles.y(jloop,iloop) > KC.L, [priorParticles.x(jloop,iloop),priorParticles.y(jloop,iloop),priorParticles.h(jloop,iloop)] = calculateWallBounce(priorParticles.x(jloop,iloop),priorParticles.y(jloop,iloop),priorParticles.h(jloop,iloop),0,[],KC.L,KC.vbar); end
    end
end

% A posteriori update/Measurement update step
beta = NaN(N,1);
beta_vals = NaN(4,1);

% calculate the distances between sensors and robots
d_1A = sqrt((priorParticles.x(1,:) -2*KC.L).^2 + priorParticles.y(1,:).^2);
d_1B = sqrt((priorParticles.x(2,:) -2*KC.L).^2 + priorParticles.y(2,:).^2);
d_2A = sqrt((priorParticles.x(1,:) -2*KC.L).^2 + (KC.L-priorParticles.y(1,:)).^2);
d_2B = sqrt((priorParticles.x(2,:) -2*KC.L).^2 + (KC.L-priorParticles.y(2,:)).^2);
d_3A = sqrt(priorParticles.x(1,:).^2 + (KC.L-priorParticles.y(1,:)).^2);
d_3B = sqrt(priorParticles.x(2,:).^2 + (KC.L-priorParticles.y(2,:)).^2);
d_4A = sqrt(priorParticles.x(1,:).^2 + (priorParticles.y(1,:)).^2);
d_4B = sqrt(priorParticles.x(2,:).^2 + (priorParticles.y(2,:)).^2);


% if we have not measurements, just resample the particles with the same
% weight
if (isinf(sens(1)) && isinf(sens(2)) && isinf(sens(3)) && isinf(sens(4)))
    beta = ones(size(beta))/N;
    
else
    % if we have measurements, then we should calculate the weight beta.
    for i = 1:N
        % For all possible cases of available sensors, generate corresponding
        % beta
        if isinf(sens(1)), beta_vals(1) = 1;
        else, beta_vals(1)  =  (1-KC.sbar)*getwProbability(sens(1)-d_1A(i), KC.wbar) + KC.sbar*getwProbability(sens(1) - d_1B(i), KC.wbar);
        end

        if isinf(sens(2)), beta_vals(2) = 1;
        else, beta_vals(2) =  (1-KC.sbar)*getwProbability(sens(2)-d_2A(i), KC.wbar) + KC.sbar*getwProbability(sens(2) - d_2B(i), KC.wbar);
        end

        if isinf(sens(3)), beta_vals(3) = 1;
        else, beta_vals(3) = (1-KC.sbar)*getwProbability(sens(3)-d_3B(i), KC.wbar) + KC.sbar*getwProbability(sens(3) - d_3A(i), KC.wbar);
        end
        if isinf(sens(4)), beta_vals(4) = 1;
        else, beta_vals(4) = (1-KC.sbar)*getwProbability(sens(4)-d_4B(i), KC.wbar) + KC.sbar*getwProbability(sens(3) - d_4A(i), KC.wbar);
        end
        beta(i) = beta_vals(1).* beta_vals(2).*beta_vals(3).*beta_vals(1);

    end
    % normalize beta
    beta= beta./sum(beta);
    
end    

% if beta have zero elments, then we get zero likelihood, we should cut off
% the zero elements in beta and the particle samples
zeor_index = find(beta==0);
priorParticles.x(:,zeor_index) = [];
priorParticles.y(:,zeor_index) = [];
priorParticles.h(:,zeor_index) = [];
beta(zeor_index) = [];

% Resample, get postParticles structure
betaSum = cumsum(beta);
r = rand(1,N);
for iloop = 1:N
    index = find(betaSum>r(iloop),1,'first');
    postParticles.x(:,iloop) = priorParticles.x(:,index);
    postParticles.y(:,iloop) = priorParticles.y(:,index);
    postParticles.h(:,iloop) = priorParticles.h(:,index);
end

% roughening
K = 0.03; %tuning parameter
d = 6;
Ei = (max([postParticles.x;postParticles.y;postParticles.h],[],2)-min([postParticles.x;postParticles.y;postParticles.h],[],2));
sigmai = K*Ei*N^(-1/d);
postParticles.x = max(min(postParticles.x+randn(size(postParticles.x)).*sigmai(1:2),2*KC.L),0);
postParticles.y = max(min(postParticles.y+randn(size(postParticles.y)).*sigmai(3:4),KC.L),0);
postParticles.h = wrapToPi(postParticles.h+randn(size(postParticles.h)).*sigmai(5:6));

end % end estimator

function [xBounce,yBounce,hBounce] = calculateWallBounce(x,y,h,directionWall,posWallx,posWally,vBar)
% function [xBounce,yBounce,hBounce] =
% calculateWallBounce(x,y,h,directionWall,posWallx,posWally,vBar) calculates the
% effect of the bounce on the wall. Assume perfect position mirroring along
% the wall. Variable directionWall is defined in radians, where the
% direction interval directionWall+[0,pi] is the room area. Leave posWallx,
% respectively posWally empty for walls parallel to the y or x axis.

directionPerpendicular = directionWall+pi/2;
% Compute bounce angle
if wrapToPi(h-directionPerpendicular) > 0
    alphaMin = -directionWall - h;
else
    alphaMin = +directionWall - h;
end
% Get noise realization
alphaPlus = alphaMin*(1+getVjRealization(vBar));
% Compute new heading
hBounce = wrapToPi(h + alphaMin + alphaPlus);
% Approximate position: mirror around wall
if ~isempty(posWallx); xBounce = posWallx - (x-posWallx); else, xBounce = x; end
if ~isempty(posWally); yBounce = posWally - (y-posWally); else, yBounce = y; end
end

function [vjReal] = getVjRealization(vBar)
% function [vjReal] = getVjRealization(vBar) will generate a realization of
% Vj, where Vj has the properties as described in Exercise #2 of the
% Programming Exercise, Recursive Estimation. vBar is a constant defined by
% the problem definition.
u = rand(1);
vjReal = nthroot(2*u-1,3)*vBar;
end

function [vsReal] = getVsRealization(vsBar)
% function [vsReal] = getVsRealization(vsBar) will generate a realization of
% Vs, where Vs has the properties as described in Exercise #2 of the
% Programming Exercise, Recursive Estimation. vsBar is a constant defined by
% the problem definition.
u = rand(1);
if u < 0.5
    vsReal = -vsBar + sqrt(2*u*vsBar^2);
else
    vsReal = vsBar - sqrt(-(u-1)*2*vsBar^2);
end
end

function [wProb] = getwProbability(wReal,wBar)
% function [wProb] = getwProbability(wReal) will generate the probability of
% the given realization wReal, where Wprob has the properties as described in Exercise #2 of the
% Programming Exercise, Recursive Estimation. wBar is a constant defined by
% the problem definition.
wProb = max(0,1/wBar - abs(wReal/wBar)*1/wBar);
end