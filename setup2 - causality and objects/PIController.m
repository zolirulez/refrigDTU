classdef PIController < handle
   properties
       ref  % Reference
       err  % Error
       int  % Integrated error
       mn   % Minimum output
       mx   % Maximum output
       T    % Sampling time
       P    % Proportional gain
       I    % Integral gain
       neg  % Invert feedback or not (-1 or 1)
   end
   methods
      function out = react(ctrl,meas)
          ctrl.err = ctrl.ref - meas;
          ctrl.int = ctrl.err*ctrl.T + ctrl.int;
          out = ctrl.neg*...
              min(max(ctrl.P*ctrl.err + ctrl.I*ctrl.int,ctrl.mn),ctrl.mx);
      end
      function initialize(ctrl,K,Ti,mn,mx,neg)
          ctrl.P = K;
          ctrl.I = ctrl.P/Ti;
          ctrl.mn = mn;
          ctrl.mx = mx;
          ctrl.neg = neg;
      end
   end
end