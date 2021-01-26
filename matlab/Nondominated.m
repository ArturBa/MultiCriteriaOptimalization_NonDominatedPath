classdef Nondominated
   properties
      T {mustBeNumeric}
      A 
      B 
      x0
   end
   methods
      function r = timeDepODE_F(obj, t, x, ut, u)
         u = interp1(ut, u, t);
         r = obj.A*x + obj.B*u;
      end
      
      function [t,x] = timeDepODE(obj, w, f)
          tspan = [0 obj.T];
          ut = linspace(0, obj.T, 25);
          u = 0.5 * sin(w*ut - f) + 0.5;

          opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
          [t,x] = ode45(@(t,x) obj.timeDepODE_F(t,x,ut,u), tspan, obj.x0, opts);
      end
      
      function j = J1(obj, x, ~)
          j = x^2 / 2 * obj.T;
      end
      
      function j = J2(obj, ~, u, w, f)
          j = 0.5*obj.T + 0.5*(-1/(w*f)*cos(f)*(cos(obj.T*u*w)+1) 
                + sin(f)*1/(w*f)*sin(obj.T*w*u));
      end
      
      function r = fitnessF(obj, WF)
          w = WF(1);
          f = WF(2);
          [~,x] = obj.timeDepODE(w, f);
          xT = x(end);
          uT = 0;5 * sin(w*obj.T - f) + 0.5;
          r = [obj.J1(xT, uT), obj.J2(xT, uT, w, f)];
      end
      
      function plotPareto(obj, paretoResp)
        figure();
        for i=1:size(paretoResp)
            w = paretoResp(i, 1);
            f = paretoResp(i, 2);
            [~, x] = obj.timeDepODE(w, f);
            xT = x(end);
            uT = 0.5 * sin(w*obj.T - f) + 0.5;
            plot(obj.J1(xT, uT),obj.J2(xT, uT, w, f),'r*')
            hold on;
        end
        xlabel('J1(x, u) = \int_0^T x_1(t) dt')
        ylabel('J2(x, u) = \int_0^T u(t) dt')
        title('Pareto Front u(t) = 0.5*sin(w*t-f) + 0.5')
        legend('Pareto front')
      end
      
      function r = paretoT(obj)
         fitG = @(WF)(obj.fitnessF(WF));
         r = gamultiobj(fitG, 2);
      end
      
      function r = searchBack(obj)
          r0 = obj.paretoT();
          r = r0;
          fitG = @(WF)(obj.fitnessF(WF));
          for t1_iter=obj.T:-0.01:0
              obj.T = t1_iter;
              r1 = gamultiobj(fitG, 2);
              
              [idx , ~] = ismember(r(:,:),r1(:,:), 'rows');
              r = r(idx, :);
              if(size(r, 1) == 1)
                  return;
              elseif (size(r, 1) == 0)
                   r = r1(1, :);
                   return
              end
          end
          r = r1(1, :);
      end
      
   end
end