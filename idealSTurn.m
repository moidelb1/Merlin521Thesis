function [s] = idealSTurn(psi0,V,phi,p,alt)
%Calculates the characteristics of a perfect S-turn assuming it is steady
%and coordinated.
%Inputs:
%   psi0 = initial aircraft heading [degrees]
%   V = aircraft speed, which remains constant [knots]
%   phi = bank angle of the turn [degrees]
%         If specified as scalar, both legs of S-turn use same bank angle,
%         and the first leg is a right turn.
%         If specified as a 2-element vector, uses different bank angles
%         for each leg. Bank angle < 0 denotes a left turn.
%   p = roll rate, assumed constant [deg/s]
%   alt = altitude of S-turn [ft]
%Outputs:
%   s = structure variable with attributes:
%       time = time since start of maneuver [seconds]
%       north = north displacement [ft]
%       east = east displacement [ft]
%       alt = altitude [ft]
%       V = aircraft speed [ft/s]
%       phi = aircraft roll/bank angle [deg]
%       psi = aircraft heading angle [deg]

if length(phi)>1
    phi1 = phi(1);
    phi2 = phi(2);
else
    phi1 = phi;
    phi2 = -phi;
end

g = 32.2;   %gravitational acceleration [ft/s^2]

%Unit Conversion
V = V*1.688;    %velocity [ft/s]

%Turn Timing
t1 = abs(phi1)/p;   %time to bank for first leg [seconds]
t2 = t1 + pi*V/g/tand(abs(phi1));   %end of first leg [seconds]
t3 = t2 + abs(phi1-phi2)/p;         %bank for second leg [seconds]
t4 = t3 + pi*V/g/tand(abs(phi2));   %end of second leg [seconds]
t5 = t4 + abs(phi2)/p;              %return to level bank [seconds]

t = 0:0.01:t5;      %time vector [seconds]

%Calculating Turn Radius
r1 = V^2/g/tand(phi1);  %ideal first leg radius [ft]
r2 = V^2/g/tand(phi2);  %ideal second leg radius [ft]

%Center Points of Turn
o1 = V*t1*[cosd(psi0),sind(psi0)] + r1*[-sind(psi0),cosd(psi0)];
o2 = o1 - V*(t3-t2)*[cosd(psi0),sind(psi0)] +...
    (r1-r2)*[-sind(psi0),cosd(psi0)];

s.time = t;                 %time vector [seconds]
s.north = zeros(size(t));   %preallocation
s.east = zeros(size(t));    %preallocation
s.alt = alt*ones(size(t));  %altitude is constant [ft/s]
s.phi = zeros(size(t));     %preallocation
s.psi = zeros(size(t));     %preallocation
s.V = V*ones(size(t));      %velocity is constant [ft/s]
for i = 1:length(t)
    if t(i) < t1    %initial banking
        s.phi(i) = sign(phi1)*t(i)*p;
        s.psi(i) = psi0;
        s.north(i) = t(i)*V*cosd(s.psi(i));
        s.east(i) = t(i)*V*sind(s.psi(i));
    elseif t(i) < t2    %first leg of turn
        s.phi(i) = phi1;
        s.psi(i) = psi0 + (t(i)-t1)*g/V*tand(s.phi(i))*180/pi;
        s.north(i) = o1(1) + r1*sind(s.psi(i));
        s.east(i) = o1(2) - r1*cosd(s.psi(i));
    elseif t(i) < t3    %bank other direction
        s.phi(i) = phi1 + sign(phi2-phi1)*(t(i)-t2)*p;
        s.psi(i) = psi0 + 180;
        s.north(i) = o1(1) - r1*sind(psi0) + (t(i)-t2)*V*cosd(s.psi(i));
        s.east(i) = o1(2) + r1*cosd(psi0) + (t(i)-t2)*V*sind(s.psi(i));
    elseif t(i) < t4    %second leg of turn
        s.phi(i) = phi2;
        s.psi(i) = psi0 + 180 + (t(i)-t3)*g/V*tand(s.phi(i))*180/pi;
        s.north(i) = o2(1) + r2*sind(s.psi(i));
        s.east(i) = o2(2) - r2*cosd(s.psi(i));
    else    %return to level flight
        s.phi(i) = phi2 - sign(phi2)*(t(i)-t4)*p;
        s.psi(i) = psi0;
        s.north(i) = o2(1) + r2*sind(psi0) + (t(i)-t4)*V*cosd(s.psi(i));
        s.east(i) = o2(2) - r2*cosd(psi0) + (t(i)-t4)*V*sind(s.psi(i));
    end
    %Correct heading to be from 0 to 360 degrees
    s.psi(i) = mod(s.psi(i),360);
end
