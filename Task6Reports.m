%% Reports

%% Part 6(a)
%
% 1.)Formulae for |AA| and |bb|
%
% $AA=\frac{1}{\Delta t}M-\theta A$
%
% $bb=(\frac{1}{\Delta t}M+\overline{\theta} A)U_n + f(\overline{\theta}t_n+\theta t_{n+1})$
%
% 2.)\ is currently used. This will result in LU factorisation and
% substitution to arive at a solution directly. The upper and lower matrices are not stored. 
% We could use an iterative
% method, by using the X0 inputted. This may result in less computations
% required to arrive at a numerical solution. Especially, if we select
% good initial conditions.
%

%% Part 6(b)
%
% 1.) 
%
%           -We have approximated a 3d problem to 2d.  
%
%           -The boundary conditions of the open ocean is an unaccurate
%           approximation, as the poluutants there will not all disperse
%           immediately. It is unreasonable hold it constant at 0.
%
%           -We have assumed that the concentration of the polutant at
%           where the pipes are located have a constant concentration. In
%           reality, the concentration at these points may not be constant,
%           as pollutants are diluted instantly when in water.
%
%           -The Neumann conditions may also not be realistic, as it may be
%           a function, and not held constant.
%
%           -The kappa term, diffusion, may not be a scalar and instead be
%           a tensor.
%
% 2.)At all three of the sample points, the concentration is below
% 1.0kgm^-1. Hence, waste concentration requirements are met. All three
% samples have concentrations that have stabilised by 20 days. Taking the
% sample at 20 days is sensible.
%
% 3.) We could compute multiple times with decreasing step size, and
% compare results. If the difference is within a decided tolerance, then we can
% terminate.
%
