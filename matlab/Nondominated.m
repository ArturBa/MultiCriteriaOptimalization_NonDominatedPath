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
      
      function [t,x] = timeDepODE(obj, f, g, h)
          tspan = [0 obj.T];
          ut = linspace(0, obj.T, 25);
          u = sign(f*ut.^2 + g*ut + h) * 0.5 + 0.5;

          opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
          [t,x] = ode45(@(t,x) obj.timeDepODE_F(t,x,ut,u), tspan, obj.x0, opts);
      end
      
      function j = J1(obj, x, ~)
          j = - x^2 / 2 * obj.T;
      end
      
      function j = J2(obj, ~, ~, f, g, h)
          r = roots([f, g, h]);
          ro = [0];
          j = 0;
          s = 0;
          if sign(h) > 0
              s = 1;
          end
          for i=1:size(r)
              if(r(i) > 0 && r(i) < obj.T)
                ro(end+1) = r(i);
              end
          end
          ro(end+1) = obj.T;
          
          for i=2:size(ro)
              j = j +  s(ro(i) - ro(i-1));
              if s == 0
                  s = 1;
              else
                  s = 0;
              end
          end
      end
      
      function r = fitnessF(obj, FGH)
          f = FGH(1);
          g = FGH(2);
          h = FGH(3);
          [~,x] = obj.timeDepODE(f, g, h);
          xT = x(end);
          uT = sign(f*obj.T^2 + g*obj.T + h)* 0.5 + 0.5;
          r = [obj.J1(xT, uT), obj.J2(xT, uT, f, g, h)];
      end
      
      function plotPareto(obj, paretoResp)
        figure();
        for i=1:size(paretoResp)
            f = paretoResp(i, 1);
            g = paretoResp(i, 2);
            h = paretoResp(i, 3);
            [~, x] = obj.timeDepODE(f, g, h);
            xT = x(end);
            plot(obj.J1(xT, 1),obj.J2(xT, 1, f, g, h),'r*')
            hold on;
        end
        xlabel('J1(x, u) = \int_0^T x_1(t) dt')
        ylabel('J2(x, u) = \int_0^T u(t) dt')
        title('Pareto Front u(t) = 0.5*sin(w*t-f) + 0.5')
        legend('Pareto front')
      end
      
      function r = paretoT(obj)
         fitG = @(FGH)(obj.fitnessF(FGH));
         r = paretosearch(fitG, 3);
      end
      
      function r = searchBack(obj)
          r0 = obj.paretoT();
          r = r0;
          fitG = @(FGH)(obj.fitnessF(FGH));
          for t1_iter=obj.T:-0.01:0
              obj.T = t1_iter;
              r1 = paretosearch(fitG, 3);
              
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