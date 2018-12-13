function sys = vfa_deriv(~,x,u) %(t,x,u)
%VFA_DERIV Returns xdot for the VFA HALE aircraft nonlinear model
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
%  thrust     = Thrust control
%  aileron_c  = Center aileron control
%  aileron_o  = Outboard aileron control
%  elevator_c = Center elevator control
%  elevator_o = Outboard elevator control
%

%% Required function defitions

rotx = @(a) [1 0 0;0 cos(a) sin(a);0 -sin(a) cos(a)];
roty = @(a) [cos(a) 0 -sin(a);0 1 0; sin(a) 0 cos(a)];
w2b = @(alpha, beta) [
    cos(alpha)*cos(beta), -cos(alpha)*sin(beta), -sin(alpha);
    sin(beta), cos(beta), 1;
    sin(alpha)*cos(beta), -sin(alpha)*sin(beta), cos(alpha)
    ];

%% Extract inputs
                  
% Start of original function
data = struct();
data.w = 300;
data.x = 30;
data.y = 20/10;
data.z = 18;

% Extract the input state parameters
V = x(1);
alpha = x(2);
h = x(3);
theta = x(4);
q = x(5);
eta = x(6);
etadot = x(7);

% Extract the input control variables
thrust = u(1);
delta2 = u(2);
delta1 = u(3);
delta3 = delta1;
deltat2 = u(4);
deltat1 = u(5);
deltat3 = deltat1;

V2dis_X = 0;
V2dis_Z = 0;
V3dis_X = 0;
V3dis_Z = 0;

% Obtain flight path angle
gamma = theta - alpha;

% Determine the air density from the standard atmosphere
[~, ~, rho, ~] = atmosphere4(h, 1);

%% Aircraft Properties

g = 32.2;

weight1wing = data.w;    % 300 lbs
mstar = weight1wing / g; % slugs
m = 3 * mstar;

ixxstar = mstar * data.x;
iyystar = mstar * data.y;
izzstar = mstar * data.z;

% Wing span and chord
s = 80;
c = 8;

areastar = s * c;

% Tail span and chord
stail = 20;
ctail = 2;

% Tail boom length
booml = 6 + 30;

% Overall dynamic pressure
rhobar = 0.5 * rho * V^2;

% Aerodynamic coefficients
kapd = 0.07;

clalpha = 2 * pi;
cld = 2;

cmdelta = -0.25;
cma = 0;
cm0 = 0.025;
cd0 = 0.007;

% Sideslip angle
beta=0;

%% Sectional velocities and aerodynamic angles

V2dis_xyz = roty(theta) * [V2dis_X, 0, V2dis_Z]';
V2dis_x = V2dis_xyz(1);
V2dis_y = V2dis_xyz(2); 
V2dis_z = V2dis_xyz(3);

V3dis_xyz = rotx(eta) * roty(theta) * [V3dis_X, 0, V3dis_Z]';
V3dis_x = V3dis_xyz(1);
V3dis_y = V3dis_xyz(2); 
V3dis_z = V3dis_xyz(3);

vcx3 = V*cos(alpha) - s/6*sin(eta)*q + V3dis_x;
vcy3 = (V*sin(alpha) + etadot*s/3*cos(eta)) * sin(eta) + V3dis_y;
vcz3 = (V*sin(alpha) + etadot*s/3*cos(eta)) * cos(eta) - etadot*s/2 + V3dis_z;

V3=sqrt(vcx3^2+vcy3^2+vcz3^2);
V1=V3;

alpha3 = atan2(vcz3, vcx3);
beta3 = asin(vcy3/V3);
alpha1 = alpha3;
beta1 = -beta3;

vbx2 = V*cos(alpha) + s/3*sin(eta)*q + V2dis_x;
vby2 = 0;
vbz2 = V*sin(alpha) + s/3*cos(eta)*etadot + V2dis_z;

V2 = sqrt(vbx2^2 + vby2^2 + vbz2^2);

alpha2 = atan2(vbz2, vbx2);
beta2 = asin(vby2/V2);

%% Determine sectional aerodynamic parameters

rhobar1 = 0.5 * rho * V1^2;
rhobar2 = 0.5 * rho * V2^2;
rhobar3 = 0.5 * rho * V3^2;

cl1 = clalpha * alpha1 + cld * delta1;
cl2 = clalpha * alpha2 + cld * delta2;
cl3 = clalpha * alpha3 + cld * delta3;

lift1 = rhobar1 *cl1 * areastar;
lift2 = rhobar2 *cl2 * areastar;
lift3 = rhobar3 *cl3 * areastar;

liftt1 = rhobar1 * clalpha * (alpha1 + deltat1) * stail * ctail;
liftt2 = rhobar2 * clalpha * (alpha2 + deltat2) * stail * ctail;
liftt3 = rhobar3 * clalpha * (alpha3 + deltat3) * stail * ctail;

drag1 = (cd0 + kapd*cl1^2) * rhobar * areastar;
drag2 = (cd0 + kapd*cl2^2) * rhobar * areastar;
drag3 = (cd0 + kapd*cl3^2) * rhobar * areastar;

%% Determine overall forces and moments

w1 = [-drag1 0 -lift1]';
w2 = [-drag2 0 -lift2]';
w3 = [-drag3 0 -lift3]';

wt1 = [0 0 -liftt1]';
wt2 = [0 0 -liftt2]';
wt3 = [0 0 -liftt3]';

p1 = rotx(eta) * w2b(alpha1, beta1) * (w1 + wt1);
p2 = w2b(alpha2, beta2) * (w2 + wt2);
p3 = rotx(-eta) * w2b(alpha3, beta3) * (w3 + wt3);

pt1 = rotx(eta) * w2b(alpha1, beta1) * wt1;
pt2 = w2b(alpha2, beta2) * wt2;
pt3 = rotx(-eta) * w2b(alpha3, beta3) * wt3;

% Overall forces

ptotal = p1+p2+p3;
wtotal = w2b(alpha, beta)' * ptotal;

lift = [0 0 -1] * wtotal;
drag = [-1 0 0] * wtotal;

% Moments

l1 = (s/2 - s/3) * sin(eta);
l2 = s/3 * sin(eta);
l3 = l1;

moment1 = rhobar1 * c * areastar * (cm0 + cmdelta*delta1 + cma*alpha1) + booml*pt1(3);
moment2 = rhobar2 * c * areastar * (cm0 + cmdelta*delta2 + cma*alpha2) + booml*pt2(3);
moment3 = rhobar3 * c * areastar * (cm0 + cmdelta*delta3 + cma*alpha3) + booml*pt3(3);

moment = moment1 + moment2 + moment3 - l1*p1(1) + l2*p2(1) - l3*p3(1);

% Dihedral angle moments

kk = 70^2;
kc = 2*sqrt(10^10)*0.707;

flapmom = -s/2*([0 0 1] * w2b(alpha3,beta3) * (w3+wt3) + mstar*g*cos(eta)*cos(theta)) - kk*eta - kc*etadot;

%% Obtain the derivative terms

%10^2 works really well 100^2 as well 10^11 for damping
Vdot = (thrust*cos(alpha) - drag)/m - g*sin(gamma);

alphadot = -(thrust*sin(alpha) + lift)/(m*V) + q + g*cos(gamma)/V;

hdot = V*sin(gamma);

thetadot = q;

c1 = 3*iyystar;
c2 = 2*izzstar - 2*iyystar + mstar*s^2/6;

qdot = (moment - 2*c2*sin(eta)*cos(eta)*etadot*q) / (c1 + c2*sin(eta)^2);

d1 = s/2*mstar * ((Vdot*sin(alpha) + V*cos(alpha)*alphadot) * cos(eta) - V*sin(alpha)*sin(eta)*etadot - 2*s/3*cos(eta)*sin(eta)*etadot^2);
d2 = (iyystar - izzstar - mstar*s^2/12) * sin(eta)*cos(eta)*q^2 - s/2*mstar*cos(eta)*V*cos(alpha)*q;
d3 = ixxstar + mstar * (s^2/4 + s^2/6*cos(eta^2));

etadotdot = (flapmom + d1 - d2) / d3;

%% Construct the system derivative
sys = [Vdot alphadot hdot thetadot qdot etadot etadotdot]';
