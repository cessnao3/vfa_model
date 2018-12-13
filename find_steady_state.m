function [ Xo, Uo, A, B ] = find_steady_state( Vo, eta, ho, fpa )
%FIND_STEADY_STATE Calculates the steady-state Xo and Uo for HALE ODE
%
% Inputs:
%  Vo  = steady-state velocity [ ft / s]
%  eta = steady-state dihedral angle [radians]
%  ho  = steady-state altitude [ft]
%  fpa = steady-state flight path angle [radians]
%
% Outputs:
%  X: Steady State
%    V     = Wind-Frame Velocity
%    alpha = Angle of Attack
%    h     = Altitude
%    theta = Pitch Angle
%    q     = Body-Frame, Longitudinal wind frame, Y rotational velocity
%    eta   = Dihedral Angle
%    etaD  = Derivative of Dihedral Angle
%
%  Uo: Control Input
%    aileron_c  = Center aileron control
%    aileron_o  = Outboard aileron control
%    elevator_c = Center elevator control
%    elevator_o = Outboard elevator control
%    thrust     = Thrust control
%
%  A: Linearized state matrix
%  B: Linearized control matrix

    % Define the initial alpha value
    alpha_guess = deg2rad(7.5) + rad2deg(eta) / 600;
    
    if nargin < 4
        fpa = 0;
    end 
    
    % Determine the input guess
    eta_deg = rad2deg(eta);
    u_guess = [0.1, deg2rad(15 + eta_deg / 3), deg2rad(5 - eta_deg / 10), -0.1, 2000]';
    
    % Set the elevator center as a known
    elev_center_guess = u_guess(3);
    
    function [X, U] = fsolve2ss( Xf )
    %FSOLVE2SS Create an internal state for fsolve to utilize
    % 
    % INPUTS
    %   X(1) = Center Aileron
    %   X(2) = Outer Aileron
    %   X(3) = Center Elevator
    %   X(4) = Outer Elevator
    %   X(5) = Thrust
    
        % Initialize variables
        X = zeros(7, 1);
        U = zeros(5, 1);
        
        % Set X
        X(1:7) = [Vo, alpha_guess, ho, fpa + alpha_guess, 0, eta, 0]';
        
        % Set U
        U(1:2) = Xf(1:2);
        U(3) = elev_center_guess;
        U(4:5) = Xf(3:4);
        
        %U = Xf;
    end

    function j = fobj( Xin )
    %FOBJ Returns the objective function (norm) for the fmincon algorithm
    %
    % INPUTS
    %   Xin(1) = Center Aileron
    %   Xin(2) = Outer Aileron
    %   Xin(3) = Outer Elevator
    %   Xin(4) = Thrust
        
        % Convert Xin to usable x, u values
        [x, u] = fsolve2ss( Xin );
        
        % Calculate xdot
        xdot = odefunc(0, x, u);
        
        % Return the norm
        j = norm(xdot);
    end

    % Create the initial guess for the standard hale_ode function
    x_start = u_guess([1:2, 4:5]);
    %x_start = u_guess;

    % Create the minimum and maximum for the xout vector
    xout_min = [-1, -1, -1, 0]';
    xout_max = [1, 1, 0, 2000]';
    
    % Run the optimization
    opts = optimoptions(@fmincon, 'Display', 'off', 'TolCon', 1e-10);
    [Xout, ~, exitflag] = fmincon(@fobj, x_start, [], [], [], [], xout_min, xout_max, [], opts);
    
    % Check for success
    if exitflag <= 0
        warning('fmincon unable to find probable solution to problem');
    end
    
    % Return the resulting values from FSOLVE
    [Xo, Uo] = fsolve2ss(Xout);
    
    % Check for linearized
    xdot_end = odefunc(0, Xo, Uo);
    if max(abs(xdot_end(1:2, 4:end))) > 1e-3
        warning('fmincon result may not be linearized %.3e', max(abs(odefunc(0, Xo, Uo))));
    end
    
    % Create the A and B linearized matrices about this design point
    threshx = ones(7, 1) * 0.01;
    threshu = ones(5, 1) * 0.01;
    
    % Calculate A
    [A, ~] = numjac(@(t, x) odefunc(t, x, Uo), 0, Xo,...
                    odefunc(0, Xo, Uo),  threshx, [], 0);

	% Calculate B
    [B, ~] = numjac(@(t, u) odefunc(t, Xo, u), 0, Uo,...
                    odefunc(0, Xo, Uo),  threshu, [], 0);
end