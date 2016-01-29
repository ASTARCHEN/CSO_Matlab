function f=benchmark_func(x,func_num)
global initial_flag
persistent fhd f_bias

if initial_flag==0
    if func_num==1      fhd=str2func('sphere_shift_func'); %[-100,100]
    elseif func_num==2  fhd=str2func('schwefel_func'); %[-100, 100]
    elseif func_num==3  fhd=str2func('rosenbrock_shift_func'); %[-100,100]
    elseif func_num==4  fhd=str2func('rastrigin_shift_func'); %[-5,5]
    elseif func_num==5  fhd=str2func('griewank_shift_func'); %[-600,600]
    elseif func_num==6  fhd=str2func('ackley_shift_func'); %[-32,32]
    elseif func_num==7  fhd=str2func('fastfractal_doubledip'); %[-1,1]    
    end
    %f_bias = [-450 -450 390 -330 -180 -140 0];
    f_bias = [0 0 0 0 0 0 0];
    %load fbias_data;
end

f=feval(fhd,x)+f_bias(func_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Unimodal%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	1.Shifted Sphere Function 
function fit=sphere_shift_func(x)
global initial_flag
persistent o
[ps,D]=size(x);
if initial_flag==0
    load sphere_shift_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-100+200*rand(1,D);
    end
    initial_flag=1;
end
x=x-repmat(o,ps,1);
fit=sum(x.^2,2);

% 2. Shifted Schwefel's Problem 2.21
function fit = schwefel_func(x)
   global initial_flag
   persistent o
   [ps, D] = size(x);
   if (initial_flag == 0)
      load schwefel_shift_func_data
      if length(o) >= D
          o=o(1:D);
      else
          o=-100+200*rand(1,D);
      end
      initial_flag = 1;
   end
   x=x-repmat(o,ps,1);
   fit = max(abs(x), [], 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Multimodal%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	3.Shifted Rosenbrock's Function
function f=rosenbrock_shift_func(x)
global initial_flag
persistent o
[ps,D]=size(x);
if initial_flag==0
    load rosenbrock_shift_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-90+180*rand(1,D);
    end
    initial_flag=1;
end
x=x-repmat(o,ps,1)+1;
f=sum(100.*(x(:,1:D-1).^2-x(:,2:D)).^2+(x(:,1:D-1)-1).^2,2);

% 4.Shifted Rastrign's Function
function f=rastrigin_shift_func(x)
global initial_flag
persistent o
[ps,D]=size(x);
if initial_flag==0
    load rastrigin_shift_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-5+10*rand(1,D);
    end
    initial_flag=1;
end
x=x-repmat(o,ps,1);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);

% 5.Shifted Griewank's Function
function f=griewank_shift_func(x)
global initial_flag
persistent o
[ps,D]=size(x);
if initial_flag==0
    load griewank_shift_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-600+1200*rand(1,D);
    end
    o=o(1:D);
    initial_flag=1;
end
x=x-repmat(o,ps,1);
f=1;
for i=1:D
    f=f.*cos(x(:,i)./sqrt(i));
end
f=sum(x.^2,2)./4000-f+1;

% 	6.Shifted Ackley's Function
function f=ackley_shift_func(x)
global initial_flag
persistent o
[ps,D]=size(x);
if initial_flag==0
    load ackley_shift_func_data
    if length(o)>=D
         o=o(1:D);
    else
         o=-30+60*rand(1,D);
    end
    initial_flag=1;
end
x=x-repmat(o,ps,1);
f=sum(x.^2,2);
f=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);

%   7. FastFractal "DoubleDip"
function f=fastfractal_doubledip(x)
global initial_flag
persistent o ff
[ps,D]=size(x);
if initial_flag==0
    load fastfractal_doubledip_data
    ff = FastFractal('DoubleDip', 3, 1, o, D);
    initial_flag=1;
end
f=ff.evaluate(x);

%%%%% end of file %%%%%

