function [max_interpolation,timedelay] = delay(X,ts)

% X is the result of the xcorr function. ts is sampling time (= 1/fs).

% This function finds the maximum of the intercorrelation function X, and
% refines the maximum by fitting a parabola passing through the maximum and the two adjacent points.
% This allows one to get a more precise time delay.

% NN is the length of the signal you want to analyze. The intercorrelation function
% has a length equal to 2*NN-1.

NN = (length(X)+1)/2;
tps = -(NN-1)*ts:ts:(NN-1)*ts;

[maxi,ind] = max(abs(X));

% parabolic interpolation around the max
x1 = ind - 1; x2 = ind; x3 = ind + 1;
y1 = X (x1);
y2 = X (x2);
y3 = X (x3);

b = ((y1 - y3)*(x1^2 - x2^2) - (y1 - y2)*(x1^2 - x3^2))/...
    ((x1 - x3)*(x1^2 - x2^2) - (x1 - x2)*(x1^2 - x3^2));
a = (y1 - y3 - b*(x1 - x3))/(x1^2 - x3^2);
c = y1 - a*x1^2 - b*x1;

ind_max_corr = -b/(2*a);
max_interpolation = a*ind_max_corr^2 + b*ind_max_corr + c;

timedelay = ts*(ind_max_corr-NN);

interpolationX = x1:1/100:x3;
interpolationY = a*interpolationX.^2 + b*interpolationX + c;

interpolationt = ts*(interpolationX-NN);

% uncomment to see the fitting result
% figure(66),
% plot(tps,X,'b');
% hold on
% % max inter
% plot(tps(ind),maxi,'*r'); % max of intercorrelation
% plot(timedelay,max_interpolation,'*g'); % refined max of intercorrelation
% plot(tps(NN),X(NN),'ok'); % interpolation
% plot(interpolationt,interpolationY,'g');
% drawnow,
% hold off
