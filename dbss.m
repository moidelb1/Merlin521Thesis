function xout = dbss(x,dband,sshape)
%Conditions the numerical data given by x using a deadband percentage and
%stick shaping percentage. For use with Merlin EIII control input data.
%
%Inputs:
%   x = numeric scalar or array
%   dband = deadband in percent of full range
%   sshape = stick shaping in percent of full range
%       Note: The input/output relationship is a blend of linear and cubic
%       functions based on this value. 100% is fully cubic and 0% is
%       completely linear.
%Outputs:
%   xout = conditioned data with the desired deadband and stick shaping
%          applied

%Create non-dimensional vector with given deadband and stick shaping
invec = -1:0.01:1;  %normalized input vector from -1 to 1
outvec = zeros(size(invec));    %preallocate output vector
for i = 1:length(invec)
    if abs(invec(i)) <= dband   %set to 0 if within deadband
        outvec(i) = 0;
    else
        sgn = sign(invec(i));   %positive or negative sign
        delt = (abs(invec(i))-dband)/(1-dband); %magnitude of input
        
        %Blend linear and cubic functions
        outvec(i) = (1-sshape)*sgn*delt + sshape*sgn*delt^3;
    end
end

%Interpolate the data
xout = interp1(invec,outvec,x);