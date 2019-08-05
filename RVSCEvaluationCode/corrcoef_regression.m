function m = corrcoef_regression(x, y, str)

%m = corrcoef_regression(x, y, str)
%
% This function compares two measurements (x and y).
% It returns a structure variable containing correlation coefficient (R)
% and linear regression fit parameters (slope a and intercept b), and
% optionally displays a plotting of x vs. y with title in str.
%
% INPUTS:
% -------
% x: reference data (manual)
% y: predicted data (automatic)
% str: string for plot title
%
% OUTPUTS:
% -------
% m: struct variable containing correlation coefficient (R) and linear regression fit parameters (slope a and intercept b).

%   Copyright: LITIS EA 4108, Université de Rouen, France
%   Author: Caroline Petitjean (caroline.petitjean@univ-rouen.fr)
%   Revision: 1.0 - Date: April 8th, 2012

%Optional plot
PLOTON = 0;

%Correlation coefficient
m.R = corr2(x,y);

%Linear regression fit
m.a = m.R*std(y)/std(x);
m.b = mean(y) - m.a*mean(x);

%Plot: Manual vs. Automatic 
if PLOTON
    figure,plot(x,y,'+');title(str);xlabel('Manual');ylabel('Automatic');
    h = (max(x)-min(x))/length(x);
    t = min(x):h:max(x);
    v = m.a*t+m.b;
    hold on,plot(t,v,'r')
end
