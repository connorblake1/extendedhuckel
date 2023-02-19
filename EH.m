classdef EH < handle
    properties
        k = 1.75; % Wolfsberg-Helmholtz constant
        slaterParams = [1.2, 2.093, .650, .975, 1.3, 1.625, 1.925, 2.275, 2.425, .733, .950]; % effective z in hydrogen wave func
        ionE=-1*[13.6 0 0;
                24.5 0 0;
                0 5.45 0;
                0 9.3 0;
                0 14 8.3;
                0 19.5 10.7;
                0 25.5 13.1;
                0 32.3 15.9];
        valenceE = -1;
        shell = [];
        funx = []; %which orbital type
        indx = []; %which atom
        zi = []; % z as a function of index
        pos = [];
        h = [];
        s = [];
        Elevels = [];
        epsSplit = .01;
        intLevels = [];
        v = [];
        occ = [];
        orbCount = -1;
    end
    methods
        function m = EH()
        end
        function molecule(m,positions, zs)
            m.zi = zs;
            m.orbCount = 0;
            m.valenceE = 0;
            for i = 1:size(m.zi,2)
                m.valenceE = m.valenceE + m.valence(m.zi(i));
                m.orbCount = m.orbCount + m.orbcount(m.zi(i));
            end
            m.shell = zeros(1,m.orbCount);
            m.funx = zeros(1,m.orbCount);
            m.indx = zeros(1,m.orbCount);
            rollOrb = 0;
            for i = 1:size(m.zi,2)
                if orbcount(m,m.zi(i)) == 1
                    m.shell(rollOrb + 1) = 1;
                    m.funx(rollOrb + 1) = 1;
                    m.indx(rollOrb + 1) = i;
                elseif orbcount(m,m.zi(i)) == 4
                    m.shell(rollOrb + 1) = 2;
                    m.shell(rollOrb + 2) = 3;
                    m.shell(rollOrb + 3) = 3;
                    m.shell(rollOrb + 4) = 3;
                    m.funx(rollOrb + 1) = 2;
                    m.funx(rollOrb + 2) = 3;
                    m.funx(rollOrb + 3) = 4;
                    m.funx(rollOrb + 4) = 5;
                    m.indx(rollOrb + 1) = i;
                    m.indx(rollOrb + 2) = i;
                    m.indx(rollOrb + 3) = i;
                    m.indx(rollOrb + 4) = i;
                end
                rollOrb = rollOrb + m.orbcount(m.zi(i));
            end
            m.pos = positions; % in angstroms
            m.pos = m.pos * 1.88973; %convert to bohr radii
            m.s = zeros(m.orbCount,m.orbCount);
        end
        function bashMatrices(m)
            for i = 1:m.orbCount % for each orbital
                for j = 1:m.orbCount
                    pi = m.pos(m.indx(i),:);
                    pj = m.pos(m.indx(j),:);
                    eval = @(x,y,z) conj(qm.CSTOgen(m.funx(i),x,y,z,pi(1),pi(2),pi(3),m.slaterParams(m.zi(m.indx(i))))).*qm.CSTOgen(m.funx(j),x,y,z,pj(1),pj(2),pj(3),m.slaterParams(m.zi(m.indx(j))));
                    m.s(i,j) = integral3(eval,-25,25,-25,25,-25,25,'RelTol',1e-15);
                    %fprintf("%d %d: %d %d %d %d %f\n",i,j,funx(i),funx(j),zi(indx(i)),zi(indx(j)),s(i,j));
                end
            end
            m.s = real(m.s);
            
            % COMPUTE H Matrix h
            m.h = zeros(m.orbCount,m.orbCount);
            for i = 1:m.orbCount
                for j = 1:m.orbCount
                    if i == j
                        m.h(i,j) = m.ionE(m.zi(m.indx(i)),m.shell(i));
                    else
                        m.h(i,j) = m.k * m.s(i,j) * .5 * (m.ionE(m.zi(m.indx(i)),m.shell(i))+m.ionE(m.zi(m.indx(j)),m.shell(j)));
                    end
                end
            end
            m.h = m.h/27.2; %in atomic units now

            % SOLVE FOR ORBITAL ENERGIES/EIGENVECTORS Elevels
            syms E;
            eqn1 = det(m.h-E*m.s) == 0;
            m.Elevels = round(double(sort(vpa(solve(eqn1,E)))),4);
            % calculates energy level degeneracy
            m.epsSplit = .01;
            m.intLevels = zeros(m.valenceE,1);
            m.intLevels(1) = 1;
            for i = 2:m.orbCount
                if m.Elevels(i) > m.Elevels(i-1) + m.epsSplit
                    m.intLevels(i) = m.intLevels(i-1) + 1;
                else
                    m.intLevels(i) = m.intLevels(i-1);
                end
            end

            % calculates MO occupancy and coefficients of AO wave functions
            m.occ = zeros(m.valenceE,1);
            m.v = zeros(m.orbCount,m.orbCount);
            j = 1;
            for i = 1:max(m.intLevels)
                orbsAtE = find(m.intLevels == i);
                for l = 1:2
                    for k1 = 1:size(orbsAtE,1)
                        vhold = null(m.h-m.s*m.Elevels(orbsAtE(k1)),.001);
                        if size(vhold,2) == 1
                            m.v(:,j) = vhold;
                        else
                            m.v(:,j) = vhold(:,k1);
                        end
                        m.occ(j) = orbsAtE(k1);
                        j = j + 1;
                    end
                end
            end

        end
        function plotElectronDensity(m,pdf)
            % PLOT SINGLE ISOSURFACE
            [x,y,z] = meshgrid([-10:0.1:10]); %#ok<NBRAK2> 
            %build wavefunction
            eval = cell(m.valenceE,1);
            for j = 1:m.valenceE
                eval{j} = @(x,y,z) 0;
                for i = 1:m.orbCount % for each orbital
                    pi = m.pos(m.indx(i),:);
                    eval{j} = @(x,y,z) eval{j}(x,y,z) + m.v(i,j)*qm.CSTOgen(m.funx(i),x,y,z,pi(1),pi(2),pi(3),m.slaterParams(m.zi(m.indx(i))));
                end
            end
            psi = @(x,y,z)0;
            for j = 1:m.valenceE
                psi = @(x,y,z) psi(x,y,z) + eval{j}(x,y,z).*eval{j}(x,y,z);
            end
            isosurface(x,y,z,feval(psi,x,y,z),pdf);
        end
        function plotDensitySlice(m,xs,ys,zs,lvls)
            hold all;
            % PLOT SINGLE ISOSURFACE
            [x,y,z] = meshgrid([-10:0.1:10]); %#ok<NBRAK2> 
            %build wavefunction
            eval = cell(m.valenceE,1);
            for j = 1:m.valenceE
                eval{j} = @(x,y,z) 0;
                for i = 1:m.orbCount % for each orbital
                    pi = m.pos(m.indx(i),:);
                    eval{j} = @(x,y,z) eval{j}(x,y,z) + m.v(i,j)*qm.CSTOgen(m.funx(i),x,y,z,pi(1),pi(2),pi(3),m.slaterParams(m.zi(m.indx(i))));
                end
            end
            psi = @(x,y,z)0;
            for j = 1:m.valenceE
                psi = @(x,y,z) psi(x,y,z) + eval{j}(x,y,z).*eval{j}(x,y,z);
            end
            contourslice(x,y,z,feval(psi,x,y,z),xs,ys,zs,lvls)
            view(3);
        end
        function luxMode(~)
            set(gcf,'Color',[0 0 0]);
            set(gca,'Color',[0 0 0]);
            set(gca,'XColor',[1 0 0]);
            set(gca,'YColor',[0 1 0]);
            set(gca,'ZColor',[0 0 1]);
            set(gca,'XLimMode','manual');
            set(gca,'YLimMode','manual');
            set(gca,'ZLimMode','manual');
            set(gca,'XLim',[-8 8]);
            set(gca,'YLim',[-8 8]);
            set(gca,'ZLim',[-8 8]);
            set(gca,'Projection','perspective');
            %then set to orbit
        end
        function num = valence(~,zin)
            if zin <= 2
                num = zin;
            elseif zin <= 10
                num = zin-2;
            elseif zin <= 18
                num = zin-10;
            end
        end
        function num = orbcount(~,zin)
            if zin <= 2
                num = 1;
            elseif zin <= 10
                num = 4;
            else
                num = 9;
            end
        end
    end
end