function [ bestever, bestSwarm ] = cso( varargin )

%  Implementation of a competitive swarm optimizer (CSO) for large scale optimization
%
%  See the details of CSO in the following paper
%  R. Cheng and Y. Jin, A Competitive Swarm Optimizer for Large Scale Optmization,
%  IEEE Transactions on Cybernetics, 2014
%
%  The test instances are the CEC'08 benchmark functions for large scale optimization
%
%  The source code CSO is implemented by Ran Cheng
%
%  If you have any questions about the code, please contact:
%  Ran Cheng at r.cheng@surrey.ac.uk
%  Prof. Yaochu Jin at yaochu.jin@surrey.ac.uk
%
%  Code update on 2018-03-18
%  A.Star, Snowland Co. Ltd , chenxiaolong12315@163.com
%
% [res, bestSwarm] = cso(fitness_function, XRRmin, XRRmax, maxfe)
% input args
%    fitness: a function handle to calculate the fitness of swarms(popsize x dim matrix)
%    popsize: integer, the popsize of swarm
%    bound_min, bound_max: the matrix(1 x dim) is the bounary of the swarms
%    maxfe: a integer, number of iterations
%    debug: default 0, display bestfitness in each iteration
% output args
% Example:
%    fun=inline('sqrt(x(:,1).^2+x(:,2).^2)','x');
%    [res, swarm] = cso(fun,100,[-8,-8],[8,8],1000,0.2);
%    xi = -8:0.05:8;
%    yi = -8:0.05:8;
%    [x,y] = meshgrid(xi, yi);
%    z = sqrt(x.^2+y.^2);
%    mesh(x,y,z);
%    hold on
%    plot3(swarm(1),swarm(2), res, 'r+')

[fitness_function,m, bound_min, bound_max, maxfe, phi, debug] = checkInput(varargin);


XRRmin = repmat(bound_min, m, 1);
XRRmax = repmat(bound_max, m, 1);
[~, d] = size(XRRmin);
%% initialization
p = XRRmin + (XRRmax - XRRmin) .* rand(m, d);
v = zeros(m,d);
bestever = 1e200;
fitness = fitness_function(p);
FES = m;
gen = 0;
bestSwarm = p(1,:);
ceil_half_m = ceil(m/2);
%% main loop
while(FES < maxfe)
    
    
    % generate random pairs
    rlist = randperm(m);
    rpairs = [rlist(1:ceil_half_m); rlist(floor(m/2) + 1:m)]';
    
    % calculate the center position
    center = ones(ceil_half_m,1)*mean(p);
    
    % do pairwise competitions
    mask = (fitness(rpairs(:,1)) > fitness(rpairs(:,2)));
    losers = mask.*rpairs(:,1) + ~mask.*rpairs(:,2);
    winners = ~mask.*rpairs(:,1) + mask.*rpairs(:,2);
    
    
    %random matrix
    randco1 = rand(ceil_half_m, d);
    randco2 = rand(ceil_half_m, d);
    randco3 = rand(ceil_half_m, d);
    
    % losers learn from winners
    v(losers,:) = randco1.*v(losers,:) ...,
        + randco2.*(p(winners,:) - p(losers,:)) ...,
        + phi*randco3.*(center - p(losers,:));
    p(losers,:) = p(losers,:) + v(losers,:);
    
    % boundary control
    for i = 1:ceil_half_m
        p(losers(i),:) = max(p(losers(i),:), XRRmin(losers(i),:));
        p(losers(i),:) = min(p(losers(i),:), XRRmax(losers(i),:));
    end
    
    
    % fitness evaluation
    fitness(losers,:) = fitness_function(p(losers,:));
    [min_fitness, min_fitness_ind] = min(fitness);
    if bestever > min_fitness
        bestever = min_fitness;
        bestSwarm = p(min_fitness_ind,:);
    end
    if debug
        fprintf('Best fitness: %e\n', bestever);
    end
    FES = FES + ceil_half_m;
    
    gen = gen + 1;
end
end


function [fitness_function,popsize, XRRmin, XRRmax, maxfe, phi, debug] = checkInput(vars)
if nargin == 7
    debug = vars{7};
else
    debug = 0;
end
fitness_function = vars{1};
popsize = vars{2};
XRRmin = vars{3};
XRRmax = vars{4};
maxfe = vars{5} * length(XRRmin);
phi = vars{6};
end