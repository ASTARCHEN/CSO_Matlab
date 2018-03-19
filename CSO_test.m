%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Implementation of a competitive swarm optimizer (CSO) for large scale optimization
%%
%%  See the details of CSO in the following paper
%%  R. Cheng and Y. Jin, A Competitive Swarm Optimizer for Large Scale Optmization,
%%  IEEE Transactions on Cybernetics, 2014
%%
%%  The test instances are the CEC'08 benchmark functions for large scale optimization
%%
%%  The source code CSO is implemented by Ran Cheng 
%%
%%  If you have any questions about the code, please contact: 
%%  Ran Cheng at r.cheng@surrey.ac.uk 
%%  Prof. Yaochu Jin at yaochu.jin@surrey.ac.uk
%% 
%%  Code update on 2018-03-18 
%%  A.Star, Snowland Co. Ltd , chenxiaolong12315@163.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath(pwd));


global initial_flag

%d: dimensionality
d = 1000;
%maxfe: maximal number of fitness evaluations
maxitor = 500;
%runnum: the number of trial runs
runnum = 1;

results = zeros(6,runnum);
 

% The frist six benchmark functions in  CEC'08  test suite.
% Function 7 is excluded because of the following error thrown by the test suite:
% 'Undefined function 'FastFractal' for input arguments of type 'char'.'
for funcid = 1 : 6
    n = d;
    initial_flag = 0;
    
    switch funcid
        case 1

            % lu: define the upper and lower bounds of the variables
            lu = [-100 * ones(1, n); 100 * ones(1, n)];

        case 2

            lu = [-100 * ones(1, n); 100 * ones(1, n)];
        case 3
            lu = [-100 * ones(1, n); 100 * ones(1, n)];
        case 4
            lu = [-5 * ones(1, n); 5 * ones(1, n)];
        case 5
            lu = [-600* ones(1, n); 600 * ones(1, n)];
        case 6
            lu = [-32 * ones(1, n); 32 * ones(1, n)];
%         case 7
%             lu = [-1 * ones(1, n); 1 * ones(1, n)];
    end
        
    %phi setting (the only parameter in CSO, please SET PROPERLY)
    if(funcid == 1 || funcid == 4 || funcid == 5 || funcid == 6)
        %for seperable functions
        if(d >= 2000)
            phi = 0.2;
        elseif(d >= 1000)
            phi = 0.15;
        elseif(d >=500)
            phi = 0.1;
        else
            phi = 0;
        end;
    else
        if(d >= 2000)
            phi = 0.2;
        elseif(d >= 1000)
            phi = 0.1;
        elseif(d >=500)
            phi = 0.05;
        else
            phi = 0;
        end;
    end;
    
    % population size setting
    if(d >= 5000)
        m = 1500;
    elseif(d >= 2000)
        m = 1000;
    elseif(d >= 1000)
        m = 500;
    elseif(d >= 100)
        m = 100;
    end;

% several runs
for run = 1 : runnum
    % initialization
    XRRmin = repmat(lu(1, :), m, 1);
    XRRmax = repmat(lu(2, :), m, 1);
    p = XRRmin + (XRRmax - XRRmin) .* rand(m, d);
    fitness_function = @(x)(benchmark_func(x, funcid));
    FES = m;
    gen = 0;
    tic;
    [bestever, swarm] = cso(fitness_function,m, lu(1, :), lu(2, :), maxitor, phi, 1);
    results(funcid, runnum) = bestever;
    fprintf('Run No.%d Done!\n', run); 
    disp(['CPU time: ',num2str(toc)]);
end;

end;


    

