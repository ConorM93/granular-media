%% Crystal Math Master 1

% doing coupled ODE similations
% hexagonal packing on hexagonal form

% indexing conventions: rows are in m, diagonal columns are n. u(m,n)
%     n-1 n n+1
%      / / /   m-1 
%     / / /     m
%    / / /     m+1
% n-1 n n+1
tic;
%% short guide:

% to change size, set L (length of side of hexagon)
% to change striker setup: change nStriker (# of strikers, line ~30) and
% striker positions sp, incident vectors sv (line ~50-60), and initial
% condition w0 right after.

%% params

% set up index matrix: stores the index in the vector w of ever sphere.
% every sphere has 4 entries, and index(m,n) has the has the w-index of the
% last of the four entries.

nu = .3; % Poisson's ratio
E = 193e9; % Young's modulus
massDens = 8*10^3;

global L num nStriker R A B Ma index sp sv scale;
scale = 1; %%DONT CHANGE THIS - I CURRENTLY DONT HAVE IT WORKING FOR NEQ 1
L = 21*scale;
nStriker = 2*scale;
R = (9.025 * 10^-3)/scale;
Ma = massDens*(4*pi*R^3/3);
num = (3*L*(L-1) + 1); tot = num + nStriker;
A = (2/3) * ( (1-nu^2)/E )^-1 * sqrt(R/2); %presimplified
B = (2/3) * ( (1-nu^2)/E )^-1 * sqrt(R);

index = zeros(2*L-1);
k = 1;
for m = 1:L
    nEnd = L+m-1;
    index(m,1:nEnd) = (k:k+nEnd-1);
    k = k+nEnd;
end
for m = L+1:2*L-1
    n1 = m - L + 1;
    index(m,n1:end) = (k:k+(2*L-1 - n1));
    k = k+(2*L-1 - n1)+1;
end


index = 4.*index;
w0 = zeros(tot, 1);


% 1) With scale 1, one strikers at 1,11:
%for i = 1:scale
%    sp(2*i-1) = 1;
%    sp(2*i) = L/2;
%    sv(2*i-1) = 0;
%    sv(2*i) = -1;
%    w0(4*(num + scale)-1) = 0;
%    w0(4*(num + scale)) = -.4*1*scale*scale;
%end %SCALE CODE NOT DONE YET
sp = [1, 11]; 
sv = [sqrt(3)/2,-1/2; -sqrt(3)/2, -1/2];


% initally at rest, striker downward velocity in y
w0(4*(num + 1)-1) = .4*sqrt(3)/2;
w0(4*(num + 1)) = -.4*1/2;
w0(4*(num + 2)-1) = -.4*sqrt(3)/2;
w0(4*(num + 2)) = -.4*1/2;


% % 2) one striker in middle of top left and top right edges:
% sp = [ceil(L/2), 1; ceil(L/2), L+ceil(L/2)-1];
% sv = [sqrt(3)/2,-1/2; -sqrt(3)/2,-1/2];
% 
% % initally at rest, striker downward velocity in y
% w0(4*(num + 1)-1) = .4*sqrt(3)/2;
% w0(4*(num + 1)) = -.4*1/2;
% w0(4*(num + 2)-1) = -.4*sqrt(3)/2;
% w0(4*(num + 2)) = -.4*1/2;


tspan = 0:0.000005:0.004; 
[t,w] = ode113(@gran_2d_hex_hex,tspan,w0);
toc
Tstep = numel(t);

%% plotting
plotIndex = zeros(2*L-1);
k = 1;
for m = 1:L
    nEnd = L+m-1;
    plotIndex(m,1:nEnd) = (k:k+nEnd-1);
    k = k+nEnd;
end
for m = L+1:2*L-1
    n1 = m - L + 1;
    nEnd = 2*L-2 - (m-L-1);
    plotIndex(m,1:nEnd) = (k:k+nEnd-1);
    k = k+nEnd;
end
%% plotting vars
v = zeros(Tstep,tot);
for i = 1:Tstep
    for j = 1:tot
        v(i,j) = norm([w(i,4*j-1),w(i,4*j)]);
    end
end
a = zeros(Tstep-1,tot);
for i = 1:Tstep-1
    for j = 1: tot
        a(i,j) = (w(i+1,4*j)-w(i,4*j))/t(2);
    end
end

%% Ignore for now: put things at right format
% VelSC = zeros(Tstep,2*L-1, 2*L-1);
% for tt = 1:Tstep
%     k=1;
%     for m = 1:L
%         nEnd = L+m-1;
%         VelSC(tt,m,1:nEnd) = v(tt,(k:k+nEnd-1));
%         k = k+nEnd;
%     end
%     for m = L+1:2*L-1
%         n1 = m - L + 1;
%         VelSC(tt,m,n1:end) = v(k:k+(2*L-1 - n1));
%         k = k+(2*L-1 - n1)+1;
%     end
% end


% VelSC is a 3D matrix, where each entry VelSC{t} is a matrix consisting of the
% velocities of all beads in the hexagonal lattice at time t.



%% Ignore for now
%
% NSc = L;
% ScNtop = 2*L-1;
% 
% Sq       = 1/2*sqrt(3);
% Xh       = [Sq 0 -Sq -Sq 0 Sq];
% Yh       = [1/2 1 1/2 -1/2 -1 -1/2];
% WaveVid2 = VideoWriter('Velocity Profile For ODE model V1.avi');
% WaveVid2.FrameRate = 5;open(WaveVid2);
% axis tight
% set(gca,'nextplot','replacechildren');
% mov(1:Tstep)=struct('cdata',[],'colormap',[]);
% 
% 
% for tt = 1:Tstep-1
%     
%     close all
%     figure(2)
%     MatSc  = squeeze(VelSC(tt, :,:)); 
%     
%     for i = 1:ScNtop
%         
%             for h = 1:NSc+1-i
%                 
%                     patch(Xh+2*Sq*h+Sq*i,Yh+3/2*(i-1),MatSc(ScNtop-(i-1),h))
%                     hold on
%                     
%             end         
%             for h = i:NSc
%                 
%                     patch(Xh+2*Sq*(h-i+1)+Sq*i,Yh-3/2*(i-1),MatSc(ScNtop+(i-1),h))
%                     hold on
%                     
%             end
%             
%     end
%     
%     axis off
%     c = colorbar;
%     c.Label.String = 'Velocity';
%     currFrame = getframe(gcf, [0 0 560 420]);
%     writeVideo(WaveVid2,currFrame);
%     
% end
% 
% close(WaveVid2)


%% evolution plot with circles

CenterPosX = zeros(size(plotIndex));
CenterPosY = zeros(size(plotIndex));
for i =1:num
    [rowPos,colPos] = find(plotIndex==i);
    CenterPosY(rowPos,colPos) = -sqrt(3)*(rowPos-1);
    if (rowPos<=L)  
        CenterPosX(rowPos,colPos) = -1*(rowPos-1)+2*(colPos-1);
    else
        CenterPosX(rowPos,colPos) = -1*(2*L-1-rowPos)+2*(colPos-1);
    end
end

%%
tPlot = (1:5:51);

figure(41);
for countPlot = 1:length(tPlot)
    %figure (countPlot+30);
    subplot(2,5,countPlot);
    hold on;
    for i = 1:num
        [rowPos,colPos] = find(plotIndex==i);
        rho = 0:0.1:1;
        theta = (0:1:360)*pi/180;
        [th,r] = meshgrid(theta,rho);
        Z = 0*th+0*r+v(tPlot(countPlot),i);
     %   Z = 0*th+0*r+material(i);
        surf(CenterPosX(rowPos,colPos)+r.*cos(th),CenterPosY(rowPos,colPos)+r.*sin(th),Z,'edgecolor','none')
        colormap parula
        %caxis([0,max(v(tPlot(countPlot),1:totConfig))]);
        %colorbar
    end
    axis off;
    str=sprintf('t = %d',tPlot(countPlot)/length(t)*max(t));
    title(str)
    hold off;
end


%%

tPlot = (11:20); tPlot = 10*tPlot + 50;

figure(42);
for countPlot = 1:length(tPlot)
    %figure (countPlot+30);
    subplot(2,5,countPlot);
    hold on;
    for i = 1:num
        [rowPos,colPos] = find(plotIndex==i);
        rho = 0:0.1:1;
        theta = (0:1:360)*pi/180;
        [th,r] = meshgrid(theta,rho);
        Z = 0*th+0*r+v(tPlot(countPlot),i);
     %   Z = 0*th+0*r+material(i);
        surf(CenterPosX(rowPos,colPos)+r.*cos(th),CenterPosY(rowPos,colPos)+r.*sin(th),Z,'edgecolor','none')
        colormap parula
        %caxis([0,max(v(tPlot(countPlot),1:totConfig))]);
        %colorbar
    end
    axis off;
    str=sprintf('t = %d',tPlot(countPlot)/length(t)*max(t));
    title(str)
    hold off;
end

%%

tPlot = (21:30); tPlot = 10*tPlot + 50;

figure(43);
for countPlot = 1:length(tPlot)
    %figure (countPlot+30);
    subplot(2,5,countPlot);
    hold on;
    for i = 1:num
        [rowPos,colPos] = find(plotIndex==i);
        rho = 0:0.1:1;
        theta = (0:1:360)*pi/180;
        [th,r] = meshgrid(theta,rho);
        Z = 0*th+0*r+v(tPlot(countPlot),i);
     %   Z = 0*th+0*r+material(i);
        surf(CenterPosX(rowPos,colPos)+r.*cos(th),CenterPosY(rowPos,colPos)+r.*sin(th),Z,'edgecolor','none')
        colormap parula
        %caxis([0,max(v(tPlot(countPlot),1:totConfig))]);
        %colorbar
    end
    axis off;
    str=sprintf('t = %d',tPlot(countPlot)/length(t)*max(t));
    title(str)
    hold off;
end

%%

tPlot = (1:10); tPlot = 5*tPlot;

figure(35);
for countPlot = 1:length(tPlot)
    %figure (countPlot+30);
    subplot(2,5,countPlot);
    hold on;
    for i = 1:num
        [rowPos,colPos] = find(plotIndex==i);
        rho = 0:0.1:1;
        theta = (0:1:360)*pi/180;
        [th,r] = meshgrid(theta,rho);
        Z = 0*th+0*r+v(tPlot(countPlot),i);
     %   Z = 0*th+0*r+material(i);
        surf(CenterPosX(rowPos,colPos)+r.*cos(th),CenterPosY(rowPos,colPos)+r.*sin(th),Z,'edgecolor','none')
        colormap parula
        %caxis([0,max(v(tPlot(countPlot),1:totConfig))]);
        %colorbar
    end
    axis off;
    str=sprintf('t = %d',tPlot(countPlot)/length(t)*max(t));
    title(str)
    hold off;
end