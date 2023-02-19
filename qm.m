classdef qm
    methods(Static)
        function psic = CSTOgen(ind,x,y,z,a,b,c,k)
            if ind == 1
                psic = qm.CSTO1s(x,y,z,a,b,c,k);
            elseif ind == 2
                psic = qm.CSTO2s(x,y,z,a,b,c,k);
            elseif ind == 3
                psic = qm.CSTO2px(x,y,z,a,b,c,k);
            elseif ind == 4
                psic = qm.CSTO2py(x,y,z,a,b,c,k);
            elseif ind == 5
                psic = qm.CSTO2pz(x,y,z,a,b,c,k);
            end
        end
        function psic = CSTO1s(x,y,z,a,b,c,k)
            r = qm.makeOffsetR(a,b,c);
            phi = qm.makeOffsetP(a,b,c);
            theta = qm.makeOffsetT(a,b,c);
            psic = qm.STO1s(r(x,y,z),phi(x,y,z),theta(x,y,z),k);
        end
        function psic = CSTO2s(x,y,z,a,b,c,k)
            r = qm.makeOffsetR(a,b,c);
            phi = qm.makeOffsetP(a,b,c);
            theta = qm.makeOffsetT(a,b,c);
            psic = qm.STO2s(r(x,y,z),phi(x,y,z),theta(x,y,z),k);
        end
        function psic = CSTO2px(x,y,z,a,b,c,k)
            r = qm.makeOffsetR(a,b,c);
            phi = qm.makeOffsetP(a,b,c);
            theta = qm.makeOffsetT(a,b,c);
            psic = qm.STO2px(r(x,y,z),phi(x,y,z),theta(x,y,z),k);
        end
        function psic = CSTO2py(x,y,z,a,b,c,k)
            r = qm.makeOffsetR(a,b,c);
            phi = qm.makeOffsetP(a,b,c);
            theta = qm.makeOffsetT(a,b,c);
            psic = qm.STO2py(r(x,y,z),phi(x,y,z),theta(x,y,z),k);
        end
        function psic = CSTO2pz(x,y,z,a,b,c,k)
            r = qm.makeOffsetR(a,b,c);
            phi = qm.makeOffsetP(a,b,c);
            theta = qm.makeOffsetT(a,b,c);
            psic = qm.STO2pz(r(x,y,z),phi(x,y,z),theta(x,y,z),k);
        end
        function psi = STO1s(r,p,t,eta) %https://www.theochem.ru.nl/~pwormer/Knowino/knowino.org/wiki/Slater_orbital.html
            psi = 1/sqrt(pi)*eta.^1.5.*exp(-r.*eta);
        end
        function psi = STO2s(r,p,t,eta)
            psi = sqrt(eta.^5/3/pi)*r.*exp(-eta*r);
        end
        function psi = STO2px(r,p,t,eta)
            psi = sqrt(eta.^5/pi)*r.*sin(p).*cos(t).*exp(-eta*r);
        end
        function psi = STO2py(r,p,t,eta)
            psi = sqrt(eta.^5/pi)*r.*sin(p).*sin(t).*exp(-eta*r);
        end
        function psi = STO2pz(r,p,t,eta)
            psi = sqrt(eta.^5/pi)*r.*cos(p).*exp(-eta*r);
        end
        
        function psic = Cggen(i,x,y,z,a,b,c,k)
            if i == 1
                psic = qm.Cg1s(x,y,z,a,b,c,k);
            elseif i == 2
                psic = qm.Cg2s(x,y,z,a,b,c,k);
            elseif i == 3
                psic = qm.Cg2px(x,y,z,a,b,c,k);
            elseif i == 4
                psic = qm.Cg2py(x,y,z,a,b,c,k);
            elseif i == 5
                psic = qm.Cg2pz(x,y,z,a,b,c,k);
            end
        end
        function psic = Cg1s(x,y,z,a,b,c,k)
            r = qm.makeOffsetR(a,b,c);
            phi = qm.makeOffsetP(a,b,c);
            theta = qm.makeOffsetT(a,b,c);
            psic = qm.oneS(r(x,y,z),phi(x,y,z),theta(x,y,z),k);
        end
        function psic = Cg2s(x,y,z,a,b,c,k)
            r = qm.makeOffsetR(a,b,c);
            phi = qm.makeOffsetP(a,b,c);
            theta = qm.makeOffsetT(a,b,c);
            psic = qm.twoS(r(x,y,z),phi(x,y,z),theta(x,y,z),k);
        end
        function psic = Cg2px(x,y,z,a,b,c,k)
            r = qm.makeOffsetR(a,b,c);
            phi = qm.makeOffsetP(a,b,c);
            theta = qm.makeOffsetT(a,b,c);
            psic = qm.twoPx(r(x,y,z),phi(x,y,z),theta(x,y,z),k);
        end
        function psic = Cg2py(x,y,z,a,b,c,k)
            r = qm.makeOffsetR(a,b,c);
            phi = qm.makeOffsetP(a,b,c);
            theta = qm.makeOffsetT(a,b,c);
            psic = qm.twoPy(r(x,y,z),phi(x,y,z),theta(x,y,z),k);
        end
        function psic = Cg2pz(x,y,z,a,b,c,k)
            r = qm.makeOffsetR(a,b,c);
            phi = qm.makeOffsetP(a,b,c);
            theta = qm.makeOffsetT(a,b,c);
            psic = qm.twoPz(r(x,y,z),phi(x,y,z),theta(x,y,z),k);
        end
        function psi = oneS(r,p,t,z)
            psi = 1/sqrt(pi)*z.^1.5.*exp(-r.*z);
        end
        function psi = twoS(r,p,t,z)
            psi = (1/sqrt(32*pi))*(z.^1.5).*(2-z.*r).*exp(-.5.*z.*r);
        end
        function psi = twoPx(r,p,t,z)
            psi = 1/sqrt(32*pi)*(z.^1.5).*z.*r.*exp(-.5*z.*r).*cos(p);
        end
        function psi = twoPy(r,p,t,z)
            psi = 1/sqrt(64*pi)*(z.^1.5).*z.*r.*exp(-.5.*z.*r).*sin(p).*exp(1i*t);
        end
        function psi = twoPz(r,p,t,z)
            psi = 1/sqrt(64*pi)*(z.^1.5).*z.*r.*exp(-.5.*z.*r).*sin(p).*exp(-1i*t);
        end
        
        function rs = makeOffsetR(a,b,c)
            rs = @r;
            function rout = r(x,y,z)
                rout = sqrt((x-a).^2 + (y-b).^2 + (z-c).^2);
            end 
        end
        function ps = makeOffsetP(a,b,c)
            ps = @phi;
            function phiout = phi(x,y,z)
                phiout = atan2(sqrt((x-a).^2+(y-b).^2),z-c);
            end
        end
        function ts = makeOffsetT(a,b,c)
            ts = @theta;
            function thetaout = theta(x,y,z)
                thetaout = atan2(y-b,x-a);
            end
        end
    end
end
