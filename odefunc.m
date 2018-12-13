% HALE Flexible Aircraft Research
% Generic ODE Function for HALE
% Ian O'Rourke

function [ xdot ] = odefunc( t, x, u )
%ODEFUNC Returns xdot for the VFA HALE aircraft, based on selected model
%
% u:
%  V     = Wind-Frame Velocity
%  alpha = Angle of Attack
%  h     = Altitude
%  theta = Pitch Angle
%  q     = Body-Frame, Longitudinal wind frame, Y rotational velocity
%  eta   = Dihedral Angle
%  etaD  = Derivative of Dihedral Angle
%
% Input u: (Controls)
%  aileron_c  = Center aileron control
%  aileron_o  = Outboard aileron control
%  elevator_c = Center elevator control
%  elevator_o = Outboard elevator control
%  thrust     = Thrust control

    % rearrange u and run
    u2 = [u(5); u(1:4)];
    xdot = vfa_deriv(t, x, u2);
end
