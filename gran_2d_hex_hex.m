function dwdt = gran_2d_hex_hex(t, w)
% given times t, position/velocity W,
% param: P(1) = L, length of one side.
%        P(2) = A, coefficient (largely constant)
%        P(3) = R, radius of spheres
%        P(4) = Ma, mass of spheres
%        P(5) = delta_0
%        P(6) = nStriker, number of strikers

%        index is a matrix holding the index of the (i,j) element in the
%        vector W / dWdt

% since the system is a second order equation, we the equation of _each_
% sphere to a 4x4 system of equations: d/dt [u, v, ud, vd].
%
%         d u_{m,n} / dt = ud_{m,n}
%         d v_{m,n} / dt = vd_{m,n}
%         d ud_{m,n} / dt = d^2 u / dt^2 = (from paper)
%         d vd_{m,n} / dt = d^2 u / dt^2 = (from paper)
% L = param(1);
% A = param(2);
% R = param(3);
% Ma = param(4);
% del = param(5);
% nStriker = param(6);

global L num nStriker R A B Ma index sp sv;

num = (3*L*(L-1) + 1); % total number of beads aka the Lth hex number + striker
tot = num + nStriker;
dwdt = zeros(4*tot,1);

%% first-order changes

% the easy entries
for i = 1:4:4*tot-3
    dwdt(i) = w(i+2);
    dwdt(i+1) = w(i+3);
end

%% spheres not on boundary

% first we do top subhalf (half, not middle row)
for m = 2:L-1
    for n = 2:L+m-2
        
        vLU = 2*R*[1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n-1)-3); w(index(m-1,n-1)-2)];
        fLU = (A/Ma)*( subplus(2*R-norm(vLU))^1.5 );
        
        vRU = 2*R*[-1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n)-3); w(index(m-1,n)-2)];
        fRU = (A/Ma)*( subplus(2*R-norm(vRU))^1.5 );
        
        vL = 2*R*[1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n-1)-3); w(index(m,n-1)-2)];
        fL = (A/Ma)*( subplus(2*R-norm(vL))^1.5 );
        
        vR = 2*R*[-1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n+1)-3); w(index(m,n+1)-2)];
        fR = (A/Ma)*( subplus(2*R-norm(vR))^1.5 );
        
        vLD = 2*R*[1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n)-3); w(index(m+1,n)-2)];
        fLD = (A/Ma)*( subplus(2*R-norm(vLD))^1.5 );
        
        vRD = 2*R*[-1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n+1)-3); w(index(m+1,n+1)-2)];
        fRD = (A/Ma)*( subplus(2*R-norm(vRD))^1.5 );
        
        % horizontal component
        dwdt(index(m,n)-1) = fLU*vLU(1)/norm(vLU) + fRU*vRU(1)/norm(vRU) + ...
            + fL*vL(1)/norm(vL) + fR*vR(1)/norm(vR) + ...
            + fLD*vLD(1)/norm(vLD) + fRD*vRD(1)/norm(vRD);
        
        % vertical component
        dwdt(index(m,n)) = fLU*vLU(2)/norm(vLU) + fRU*vRU(2)/norm(vRU) + ...
            + fL*vL(2)/norm(vL) + fR*vR(2)/norm(vR) + ...
            + fLD*vLD(2)/norm(vLD) + fRD*vRD(2)/norm(vRD);
    end
end

m = L;
for n = 2:2*L-2
    
    vLU = 2*R*[1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n-1)-3); w(index(m-1,n-1)-2)];
    fLU = (A/Ma)*( subplus(2*R-norm(vLU))^1.5 );
    
    vRU = 2*R*[-1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n)-3); w(index(m-1,n)-2)];
    fRU = (A/Ma)*( subplus(2*R-norm(vRU))^1.5 );
    
    vL = 2*R*[1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n-1)-3); w(index(m,n-1)-2)];
    fL = (A/Ma)*( subplus(2*R-norm(vL))^1.5 );
    
    vR = 2*R*[-1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n+1)-3); w(index(m,n+1)-2)];
    fR = (A/Ma)*( subplus(2*R-norm(vR))^1.5 );
    
    vLD = 2*R*[1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n)-3); w(index(m+1,n)-2)];
    fLD = (A/Ma)*( subplus(2*R-norm(vLD))^1.5 );
    
    vRD = 2*R*[-1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n+1)-3); w(index(m+1,n+1)-2)];
    fRD = (A/Ma)*( subplus(2*R-norm(vRD))^1.5 );
    
    % horizontal component
    dwdt(index(m,n)-1) = fLU*vLU(1)/norm(vLU) + fRU*vRU(1)/norm(vRU) + ...
        + fL*vL(1)/norm(vL) + fR*vR(1)/norm(vR) + ...
        + fLD*vLD(1)/norm(vLD) + fRD*vRD(1)/norm(vRD);
    
    % vertical component
    dwdt(index(m,n)) = fLU*vLU(2)/norm(vLU) + fRU*vRU(2)/norm(vRU) + ...
        + fL*vL(2)/norm(vL) + fR*vR(2)/norm(vR) + ...
        + fLD*vLD(2)/norm(vLD) + fRD*vRD(2)/norm(vRD);
    
end

% bottom subhalf
for m = L+1 : 2*L-2
    for n = m - L + 2 : 2*L-2
        a = 1; 
        
        vLU = 2*R*[1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n-1)-3); w(index(m-1,n-1)-2)];
        fLU = (A/Ma)*( subplus(2*R-norm(vLU))^1.5 );
        
        vRU = 2*R*[-1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n)-3); w(index(m-1,n)-2)];
        fRU = (A/Ma)*( subplus(2*R-norm(vRU))^1.5 );
        
        vL = 2*R*[1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n-1)-3); w(index(m,n-1)-2)];
        fL = (A/Ma)*( subplus(2*R-norm(vL))^1.5 );
        
        vR = 2*R*[-1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n+1)-3); w(index(m,n+1)-2)];
        fR = (A/Ma)*( subplus(2*R-norm(vR))^1.5 );
        
        vLD = 2*R*[1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n)-3); w(index(m+1,n)-2)];
        fLD = (A/Ma)*( subplus(2*R-norm(vLD))^1.5 );
        
        vRD = 2*R*[-1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n+1)-3); w(index(m+1,n+1)-2)];
        fRD = (A/Ma)*( subplus(2*R-norm(vRD))^1.5 );
        
        % horizontal component
        dwdt(index(m,n)-1) = fLU*vLU(1)/norm(vLU) + fRU*vRU(1)/norm(vRU) + ...
            + fL*vL(1)/norm(vL) + fR*vR(1)/norm(vR) + ...
            + fLD*vLD(1)/norm(vLD) + fRD*vRD(1)/norm(vRD);
        
        % vertical component
        dwdt(index(m,n)) = fLU*vLU(2)/norm(vLU) + fRU*vRU(2)/norm(vRU) + ...
            + fL*vL(2)/norm(vL) + fR*vR(2)/norm(vR) + ...
            + fLD*vLD(2)/norm(vLD) + fRD*vRD(2)/norm(vRD);
    end
end

%% straight-line boundary conditions (not corners)

% first row
m = 1;
for n = 2:L-1
    
    vWall = [0;-1];
    dWall = -vWall' * [w(index(m,n)-3);w(index(m,n)-2)];
    fWall = (B/Ma)*subplus(dWall)^1.5;
    
    vL = 2*R*[1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n-1)-3); w(index(m,n-1)-2)];
    fL = (A/Ma)*( subplus(2*R-norm(vL))^1.5 );
    
    vR = 2*R*[-1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n+1)-3); w(index(m,n+1)-2)];
    fR = (A/Ma)*( subplus(2*R-norm(vR))^1.5 );
    
    vLD = 2*R*[1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n)-3); w(index(m+1,n)-2)];
    fLD = (A/Ma)*( subplus(2*R-norm(vLD))^1.5 );
    
    vRD = 2*R*[-1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n+1)-3); w(index(m+1,n+1)-2)];
    fRD = (A/Ma)*( subplus(2*R-norm(vRD))^1.5 );
    
    % horizontal component
    dwdt(index(m,n)-1) = fWall*vWall(1) + ...
        + fL*vL(1)/norm(vL) + fR*vR(1)/norm(vR) + ...
        + fLD*vLD(1)/norm(vLD) + fRD*vRD(1)/norm(vRD);
    
    % vertical component
    dwdt(index(m,n)) = fWall*vWall(2) + ...
        + fL*vL(2)/norm(vL) + fR*vR(2)/norm(vR) + ...
        + fLD*vLD(2)/norm(vLD) + fRD*vRD(2)/norm(vRD);
    
end

% last row
m = 2*L-1;
for n = L+1:2*L-2
    vLU = 2*R*[1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n-1)-3); w(index(m-1,n-1)-2)];
    fLU = (A/Ma)*( subplus(2*R-norm(vLU))^1.5 );
    
    vRU = 2*R*[-1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n)-3); w(index(m-1,n)-2)];
    fRU = (A/Ma)*( subplus(2*R-norm(vRU))^1.5 );
    
    vL = 2*R*[1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n-1)-3); w(index(m,n-1)-2)];
    fL = (A/Ma)*( subplus(2*R-norm(vL))^1.5 );
    
    vR = 2*R*[-1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n+1)-3); w(index(m,n+1)-2)];
    fR = (A/Ma)*( subplus(2*R-norm(vR))^1.5 );
    
    vWall = [0;1];
    dWall = -vWall' * [w(index(m,n)-3);w(index(m,n)-2)];
    fWall = (B/Ma)*subplus(dWall)^1.5;
    
    
    % horizontal component
    dwdt(index(m,n)-1) = fLU*vLU(1)/norm(vLU) + fRU*vRU(1)/norm(vRU) + ...
        + fL*vL(1)/norm(vL) + fR*vR(1)/norm(vR) + ...
        + fWall*vWall(1);
    
    % vertical component
    dwdt(index(m,n)) = fLU*vLU(2)/norm(vLU) + fRU*vRU(2)/norm(vRU) + ...
        + fL*vL(2)/norm(vL) + fR*vR(2)/norm(vR) + ...
        + fWall*vWall(2);
    
end

% left upper boundary
n = 1;
for m = 2:L-1
    
    vWall = [sqrt(3)/2;-1/2];
    dWall = -vWall' * [w(index(m,n)-3);w(index(m,n)-2)];
    fWall = (B/Ma)*subplus(dWall)^1.5;
    
    vRU = 2*R*[-1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n)-3); w(index(m-1,n)-2)];
    fRU = (A/Ma)*( subplus(2*R-norm(vRU))^1.5 );
    
    vR = 2*R*[-1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n+1)-3); w(index(m,n+1)-2)];
    fR = (A/Ma)*( subplus(2*R-norm(vR))^1.5 );
    
    vLD = 2*R*[1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n)-3); w(index(m+1,n)-2)];
    fLD = (A/Ma)*( subplus(2*R-norm(vLD))^1.5 );
    
    vRD = 2*R*[-1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n+1)-3); w(index(m+1,n+1)-2)];
    fRD = (A/Ma)*( subplus(2*R-norm(vRD))^1.5 );
    
    % horizontal component
    dwdt(index(m,n)-1) = fWall*vWall(1) + fRU*vRU(1)/norm(vRU) + ...
        + fR*vR(1)/norm(vR) + ...
        + fLD*vLD(1)/norm(vLD) + fRD*vRD(1)/norm(vRD);
    
    % vertical component
    dwdt(index(m,n)) = fWall*vWall(2) + fRU*vRU(2)/norm(vRU) + ...
        + fR*vR(2)/norm(vR) + ...
        + fLD*vLD(2)/norm(vLD) + fRD*vRD(2)/norm(vRD);
end

% left lower boundary
% black magic
for n = 2:L-1
    m = L+n-1;
    
    vLU = 2*R*[1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n-1)-3); w(index(m-1,n-1)-2)];
    fLU = (A/Ma)*( subplus(2*R-norm(vLU))^1.5 );
    
    vRU = 2*R*[-1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n)-3); w(index(m-1,n)-2)];
    fRU = (A/Ma)*( subplus(2*R-norm(vRU))^1.5 );
    
    vR = 2*R*[-1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n+1)-3); w(index(m,n+1)-2)];
    fR = (A/Ma)*( subplus(2*R-norm(vR))^1.5 );
    
    vWall = [sqrt(3)/2;1/2];
    dWall = -vWall' * [w(index(m,n)-3);w(index(m,n)-2)];
    fWall = (B/Ma)*subplus(dWall)^1.5;
    
    vRD = 2*R*[-1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n+1)-3); w(index(m+1,n+1)-2)];
    fRD = (A/Ma)*( subplus(2*R-norm(vRD))^1.5 );
    
    % horizontal component
    dwdt(index(m,n)-1) = fLU*vLU(1)/norm(vLU) + fRU*vRU(1)/norm(vRU) + ...
        + fR*vR(1)/norm(vR) + ...
        + fWall*vWall(1) + fRD*vRD(1)/norm(vRD);
    
    % vertical component
    dwdt(index(m,n)) = fLU*vLU(2)/norm(vLU) + fRU*vRU(2)/norm(vRU) + ...
        + fR*vR(2)/norm(vR) + ...
        + fWall*vWall(2) + fRD*vRD(2)/norm(vRD);
    
end

% right upper boundary
for m = 2:L-1
    n = L+m-1;
    vLU = 2*R*[1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n-1)-3); w(index(m-1,n-1)-2)];
    fLU = (A/Ma)*( subplus(2*R-norm(vLU))^1.5 );
    
    vWall = [-sqrt(3)/2;-1/2];
    dWall = -vWall' * [w(index(m,n)-3);w(index(m,n)-2)];
    fWall = (B/Ma)*subplus(dWall)^1.5;
    
    vL = 2*R*[1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n-1)-3); w(index(m,n-1)-2)];
    fL = (A/Ma)*( subplus(2*R-norm(vL))^1.5 );
    
    vLD = 2*R*[1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n)-3); w(index(m+1,n)-2)];
    fLD = (A/Ma)*( subplus(2*R-norm(vLD))^1.5 );
    
    vRD = 2*R*[-1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n+1)-3); w(index(m+1,n+1)-2)];
    fRD = (A/Ma)*( subplus(2*R-norm(vRD))^1.5 );
    
    % horizontal component
    dwdt(index(m,n)-1) = fLU*vLU(1)/norm(vLU) + fWall*vWall(1) + ...
        + fL*vL(1)/norm(vL) + ...
        + fLD*vLD(1)/norm(vLD) + fRD*vRD(1)/norm(vRD);
    
    % vertical component
    dwdt(index(m,n)) = fLU*vLU(2)/norm(vLU) + fWall*vWall(2) + ...
        + fL*vL(2)/norm(vL) + ...
        + fLD*vLD(2)/norm(vLD) + fRD*vRD(2)/norm(vRD);
end

% right lower boundary
n = 2*L-1;
for m = L+1:2*L-2
        vLU = 2*R*[1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n-1)-3); w(index(m-1,n-1)-2)];
        fLU = (A/Ma)*( subplus(2*R-norm(vLU))^1.5 );
        
        vRU = 2*R*[-1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n)-3); w(index(m-1,n)-2)];
        fRU = (A/Ma)*( subplus(2*R-norm(vRU))^1.5 );
        
        vL = 2*R*[1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n-1)-3); w(index(m,n-1)-2)];
        fL = (A/Ma)*( subplus(2*R-norm(vL))^1.5 );
        
        vLD = 2*R*[1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n)-3); w(index(m+1,n)-2)];
        fLD = (A/Ma)*( subplus(2*R-norm(vLD))^1.5 );
        
        vWall = [-sqrt(3)/2; 1/2];
        dWall = -vWall' * [w(index(m,n)-3);w(index(m,n)-2)];
        fWall = (B/Ma)*subplus(dWall)^1.5;
        
        % horizontal component
        dwdt(index(m,n)-1) = fLU*vLU(1)/norm(vLU) + fRU*vRU(1)/norm(vRU) + ...
            + fL*vL(1)/norm(vL) + ...
            + fLD*vLD(1)/norm(vLD) + fWall*vWall(1);
        
        % vertical component
        dwdt(index(m,n)) = fLU*vLU(2)/norm(vLU) + fRU*vRU(2)/norm(vRU) + ...
            + fL*vL(2)/norm(vL) + ...
            + fLD*vLD(2)/norm(vLD) + fWall*vWall(2);
end

%% corners

% top left
m = 1; n = 1;
    vWallU = [0;-1]; vWallL = [sqrt(3)/2;-1/2];
    dWallU = -vWallU' * [w(index(m,n)-3);w(index(m,n)-2)];
    dWallL = -vWallL' * [w(index(m,n)-3);w(index(m,n)-2)];
    fWall = (B/Ma)*subplus(dWallU)^1.5 * vWallU + (B/Ma)*subplus(dWallL)^1.5 * vWallL;

        vR = 2*R*[-1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n+1)-3); w(index(m,n+1)-2)];
        fR = (A/Ma)*( subplus(2*R-norm(vR))^1.5 );
        
        vLD = 2*R*[1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n)-3); w(index(m+1,n)-2)];
        fLD = (A/Ma)*( subplus(2*R-norm(vLD))^1.5 );
        
        vRD = 2*R*[-1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n+1)-3); w(index(m+1,n+1)-2)];
        fRD = (A/Ma)*( subplus(2*R-norm(vRD))^1.5 );
        
        % horizontal component
        dwdt(index(m,n)-1) = fWall(1) + fR*vR(1)/norm(vR) + ...
            + fLD*vLD(1)/norm(vLD) + fRD*vRD(1)/norm(vRD);
        
        % vertical component
        dwdt(index(m,n)) = fWall(2) + fR*vR(2)/norm(vR) + ...
            + fLD*vLD(2)/norm(vLD) + fRD*vRD(2)/norm(vRD);

% top right
m = 1; n = L;
        vWallU = [0;-1]; vWallR = [-sqrt(3)/2;-1/2];
        dWallU = -vWallU' * [w(index(m,n)-3);w(index(m,n)-2)];
        dWallR = -vWallR' * [w(index(m,n)-3);w(index(m,n)-2)];
        fWall = (B/Ma)*subplus(dWallU)^1.5 * vWallU + (B/Ma)*subplus(dWallR)^1.5 * vWallR;
    
        vL = 2*R*[1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n-1)-3); w(index(m,n-1)-2)];
        fL = (A/Ma)*( subplus(2*R-norm(vL))^1.5 );
        
        vLD = 2*R*[1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n)-3); w(index(m+1,n)-2)];
        fLD = (A/Ma)*( subplus(2*R-norm(vLD))^1.5 );
        
        vRD = 2*R*[-1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n+1)-3); w(index(m+1,n+1)-2)];
        fRD = (A/Ma)*( subplus(2*R-norm(vRD))^1.5 );
        
        % horizontal component
        dwdt(index(m,n)-1) = fWall(1) + fL*vL(1)/norm(vR) + ...
            + fLD*vLD(1)/norm(vLD) + fRD*vRD(1)/norm(vRD);
        
        % vertical component
        dwdt(index(m,n)) = fWall(2) + fL*vL(2)/norm(vR) + ...
            + fLD*vLD(2)/norm(vLD) + fRD*vRD(2)/norm(vRD);

% middle left
m = L; n = 1;
        vWallU = [sqrt(3)/2;-1/2]; vWallD = [sqrt(3)/2;1/2];
        dWallU = -vWallU' * [w(index(m,n)-3);w(index(m,n)-2)];
        dWallD = -vWallD' * [w(index(m,n)-3);w(index(m,n)-2)];
        fWall = (B/Ma)*subplus(dWallU)^1.5 * vWallU + (B/Ma)*subplus(dWallD)^1.5 * vWallD;
        
        vRU = 2*R*[-1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n)-3); w(index(m-1,n)-2)];
        fRU = (A/Ma)*( subplus(2*R-norm(vRU))^1.5 );
        
        vR = 2*R*[-1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n+1)-3); w(index(m,n+1)-2)];
        fR = (A/Ma)*( subplus(2*R-norm(vR))^1.5 );
        
        vRD = 2*R*[-1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n+1)-3); w(index(m+1,n+1)-2)];
        fRD = (A/Ma)*( subplus(2*R-norm(vRD))^1.5 );
        
        % horizontal component
        dwdt(index(m,n)-1) = fRU*vRU(1)/norm(vRU) + fR*vR(1)/norm(vR) + ...
            + fRD*vRD(1)/norm(vRD) + fWall(1);
        
        % vertical component
        dwdt(index(m,n)) =  fRU*vRU(2)/norm(vRU) + fR*vR(2)/norm(vR) + ...
            + fRD*vRD(2)/norm(vRD) + fWall(2);

% middle right
m = L; n = 2*L-1;
        vWallU = [-sqrt(3)/2;-1/2]; vWallD = [-sqrt(3)/2;1/2];
        dWallU = -vWallU' * [w(index(m,n)-3);w(index(m,n)-2)];
        dWallD = -vWallD' * [w(index(m,n)-3);w(index(m,n)-2)];
        fWall = (B/Ma)*subplus(dWallU)^1.5 * vWallU + (B/Ma)*subplus(dWallD)^1.5 * vWallD;
        
        vLU = 2*R*[1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n-1)-3); w(index(m-1,n-1)-2)];
        fLU = (A/Ma)*( subplus(2*R-norm(vLU))^1.5 );
        
        vL = 2*R*[1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n-1)-3); w(index(m,n-1)-2)];
        fL = (A/Ma)*( subplus(2*R-norm(vL))^1.5 );
        
        vLD = 2*R*[1/2; sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m+1,n)-3); w(index(m+1,n)-2)];
        fLD = (A/Ma)*( subplus(2*R-norm(vLD))^1.5 );
        
        % horizontal component
        dwdt(index(m,n)-1) = fLU*vLU(1)/norm(vLU) + fL*vL(1)/norm(vL) + ...
            + fLD*vLD(1)/norm(vLD) + fWall(1);
        
        % vertical component
        dwdt(index(m,n)) = fLU*vLU(2)/norm(vLU) + fL*vL(2)/norm(vL) + ...
            + fLD*vLD(2)/norm(vLD) + fWall(2);

% bottom left
m = 2*L-1; n = L;
        vWallL = [sqrt(3)/2;1/2]; vWallD = [0;1];
        dWallL = -vWallL' * [w(index(m,n)-3);w(index(m,n)-2)];
        dWallD = -vWallD' * [w(index(m,n)-3);w(index(m,n)-2)];
        fWall = (B/Ma)*subplus(dWallL)^1.5 * vWallL + (B/Ma)*subplus(dWallD)^1.5 * vWallD;
        
        vLU = 2*R*[1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n-1)-3); w(index(m-1,n-1)-2)];
        fLU = (A/Ma)*( subplus(2*R-norm(vLU))^1.5 );
        
        vRU = 2*R*[-1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n)-3); w(index(m-1,n)-2)];
        fRU = (A/Ma)*( subplus(2*R-norm(vRU))^1.5 );
        
        vR = 2*R*[-1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n+1)-3); w(index(m,n+1)-2)];
        fR = (A/Ma)*( subplus(2*R-norm(vR))^1.5 );
        
        % horizontal component
        dwdt(index(m,n)-1) = fLU*vLU(1)/norm(vLU) + fRU*vRU(1)/norm(vRU) + ...
            + fR*vR(1)/norm(vR) + fWall(1) ;
        
        % vertical component
        dwdt(index(m,n)) = fLU*vLU(2)/norm(vLU) + fRU*vRU(2)/norm(vRU) + ...
            + fR*vR(2)/norm(vR) + fWall(2);

% bottom right
m = 2*L-1; n = 2*L-1;
        vWallR = [-sqrt(3)/2;1/2]; vWallD = [0;1];
        dWallR = -vWallR' * [w(index(m,n)-3);w(index(m,n)-2)];
        dWallD = -vWallD' * [w(index(m,n)-3);w(index(m,n)-2)];
        fWall = (B/Ma)*subplus(dWallR)^1.5 * vWallR + (B/Ma)*subplus(dWallD)^1.5 * vWallD;
        
        vLU = 2*R*[1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n-1)-3); w(index(m-1,n-1)-2)];
        fLU = (A/Ma)*( subplus(2*R-norm(vLU))^1.5 );
        
        vRU = 2*R*[-1/2; -sqrt(3)/2] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m-1,n)-3); w(index(m-1,n)-2)];
        fRU = (A/Ma)*( subplus(2*R-norm(vRU))^1.5 );
        
        vL = 2*R*[1; 0] + [w(index(m,n)-3); w(index(m,n)-2)] - [w(index(m,n-1)-3); w(index(m,n-1)-2)];
        fL = (A/Ma)*( subplus(2*R-norm(vL))^1.5 );
        
        % horizontal component
        dwdt(index(m,n)-1) = fLU*vLU(1)/norm(vLU) + fRU*vRU(1)/norm(vRU) + ...
            + fL*vL(1)/norm(vL) + fWall(1);
        
        % vertical component
        dwdt(index(m,n)) = fLU*vLU(2)/norm(vLU) + fRU*vRU(2)/norm(vRU) + ...
            + fL*vL(2)/norm(vL) + fWall(2);

%% strikers
% nStriker = # of strikers
% both matrices [nStriker, 2], sp has striker position, sv incident angle;

for i = 1:nStriker
    m = sp(i,1); n = sp(i,2);
    vStrike = - 2*R*sv(i,:)' + [w(4*(num+i)-3);w(4*(num+i)-2)] - [w(index(m,n)-3);w(index(m,n)-2)];
    
    dwdt(4*(num+i)-1) = (A/Ma)*( subplus(2*R-norm(vStrike))^1.5 * vStrike(1)/norm(vStrike));
    dwdt(4*(num+i)) = (A/Ma)*( subplus(2*R-norm(vStrike))^1.5 * vStrike(2)/norm(vStrike));
    
    dwdt(index(m,n)-1) = dwdt(index(m,n)-1) - (A/Ma)*( subplus(2*R - norm(vStrike))^1.5 * vStrike(1)/norm(vStrike) );
    dwdt(index(m,n)) = dwdt(index(m,n)) - (A/Ma)*( subplus(2*R - norm(vStrike))^1.5 * vStrike(2)/norm(vStrike) );
    
end


end