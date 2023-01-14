% Ivan Rendo Barreiro - QTE - R&D dopoly (Ericson and Pakes)
% Modication of Ben Moll's HJB_no_uncertainty_implicit.m code.
% with some lines from 'HJB_diffusion_implicit.m' also fron B. Moll.

tic; 
clear all; clc;

rho=0.05; 
a = 0.03;  % probability parameter
D = 10; % market size
m=20; % cost Mx^2 
I=20;

cmin = 0;
cmax = 10;
c1 = linspace(cmin,cmax,I)';
c2 = linspace(cmin,cmax,I);

cc1 = c1*ones(1,I);
cc2 = ones(I,1)*c2;

dc = (cmax-cmin)/(I-1);

maxit=1000;
crit = 10^(-6);

Delta = 1000; % mas grande, mas lento pero mas estable

dV1f = zeros(I,I);
dV2f = zeros(I,I);
dV1b = zeros(I,I);
dV2b = zeros(I,I);

%INITIAL GUESSES
profit1 = zeros(I,I);
profit2 = zeros(I,I);
for i=1:I
    for j=1:I
        profit1(i,j) = max(0,(1/9)*(D-2*c1(i)+c2(j))*(D+c1(i)+c2(j)));
        profit2(i,j) = max(0,(1/9)*(D-2*c2(j)+c1(i))*(D+c2(j)+c1(i)));
    end
end

v1 = profit1 / rho; % net present value of profits in the SS, assuming q1=q2
v2 = profit2 / rho;


% policy function
x1 = dV1f;
x2 = dV2f;

for n=1:maxit
    V1 = v1;
    V2 = v2;
    % forward difference
    dV1f(1:I-1,:) = (V1(2:I,:)-V1(1:I-1,:))/dc;
    dV2f(1:I-1,:) = (V2(2:I,:)-V2(1:I-1,:))/dc;
    
    dV1f(I,:) = dV1b(I-1,:); % (p(I)-c-0)/s ; % using FOC and q>=0 % CAMBIAAAAAAAAAAAAR
    dV2f(I,:) = dV2b(I-1,:); % (p(I)-c-0)/s; % using FOC and q>=0  % CAMBIAAAAAAAAAAAAR
    % backward difference
    dV1b(2:I,:) = (V1(2:I,:)-V1(1:I-1,:))/dc;
    dV2b(2:I,:) = (V2(2:I,:)-V2(1:I-1,:))/dc;
    
    dV1b(1,:) = dV1b(2,:); % using FOC and q>=0    % CAMBIAAAAAAAAAAAAR
    dV2b(1,:) = dV2b(2,:); % using FOC and q>=0    % CAMBIAAAAAAAAAAAAR
    
    
    % quantities and price change with forward difference
    x1f = max(a.*dV1f./m, 0);
    x2f = max(a.*dV2f./m, 0); 
    p1m = x1f;
    p1l = 1-x1f;
    p2m = x2f;
    p2l = 1-x2f;
    
    c1f = p1m.*p2m + p1m.*p2l - p1l.*p2m - p1l.*p2l;
    c2f = p2m.*p1m + p2m.*p1l - p2l.*p1m - p2l.*p1l;
    
    % quantities and price change backward difference
    x1b = max(a.*dV1b./m, 0);
    x2b = max(a.*dV2b./m, 0); 
    p1m = x1b;
    p1l = 1-x1b;
    p2m = x2b;
    p2l = 1-x2b;
    
    c1b = p1m.*p2m + p1m.*p2l - p1l.*p2m - p1l.*p2l;
    c2b = p2m.*p1m + p2m.*p1l - p2l.*p1m - p2l.*p1l; % by definition of c'
    
    %consumption and derivative of value function at steady state
    con = 0;
    dV10 = 0.5*mean(dV1b) + 0.5*mean(dV1f); % I do not have the actual value
    dV20 = 0.5*mean(dV2b) + 0.5*mean(dV2f); 
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    I1f = c1f > 0; %positive drift --> forward difference
    I2f = c2f > 0;
    I1b = c1b < 0; %negative drift --> backward difference
    I2b = c2b < 0;
    
    Iff = I1f*I2f;
    Ifb = I1f*I2b;
    Ibf = I1b*I2f;
    Ibb = I1b*I2b;
    
    I0 = (1-Iff-Ifb-Ibf-Ibb); 
    
    dV1_Upwind = dV1f.*(Iff+Ifb) + dV1b.*(Ibf+Ibb) + dV10.*I0;
    dV2_Upwind = dV2f.*(Iff+Ibf) + dV2b.*(Ifb+Ibb) + dV20.*I0;
    
    x1 = max(a.*dV1_Upwind./m, 0);
    x2 = max(a.*dV2_Upwind./m, 0);
    
    cero = zeros(I,I);
    %pi1 = (1/9).*((D-2.*cc1+cc2).*(D+cc1+cc2))-m*x1.^2;
    %pi2 = (1/9).*((D-2.*cc2+cc1).*(D+cc1+cc2))-m*x2.^2;
    bpi1 = (1/9).*((D-2.*cc1+cc2).*(D+cc1+cc2));
    bpi2 = (1/9).*((D-2.*cc2+cc1).*(D+cc1+cc2));
    
    monopoly1 = bpi2<0;
    monopoly2 = bpi1<0;
    
    pi1mono = ((D+cc1).*(D-cc1)) ./ 4;
    pi2mono = ((D+cc2).*(D-cc2)) ./ 4;
    
    bpi1(monopoly1==1) = pi1mono(monopoly1==1);
    bpi2(monopoly1==1) = cero(monopoly1==1);
    
    bpi2(monopoly2==1) = pi2mono(monopoly2==1);
    bpi1(monopoly2==1) = cero(monopoly2==1);
    
    pi1 = bpi1 - m*x1.^2;
    pi2 = bpi2 - m*x2.^2;

    
    %CONSTRUCT MATRIX A
    
    X = -min(c1b,0)/dc;
    Y = -max(c1f,0)/dc + min(c1b,0)/dc;
    Z =  max(c1f,0)/dc;
    
    updiag=0; %This is needed because of the peculiarity of spdiags.
    for j=1:I
        updiag=[updiag;Z(j,1:I-1)';0];
    end
    
    centdiag=reshape(Y',I*I,1);
    
    lowdiag=X(2:I,1);
    for j=2:I
        lowdiag=[lowdiag;0;X(j,2:I)'];
    end
    
    AA1=spdiags(centdiag,0,I*I,I*I)+spdiags([updiag;0],1,I*I,I*I)+spdiags([lowdiag;0],-1,I*I,I*I);
    
    %CONSTRUCT MATRIX A
    X = -min(c2b,0)/dc;
    Y = -max(c2f,0)/dc + min(c2b,0)/dc;
    Z =  max(c2f,0)/dc;
    
    updiag=0; %This is needed because of the peculiarity of spdiags.
    for j=1:I
        updiag=[updiag;Z(1:I-1,j);0];
    end
    
    centdiag=reshape(Y,I*I,1);
    
    lowdiag=X(2:I,1);
    for j=2:I
        lowdiag=[lowdiag;0;X(2:I,j)];
    end
    
    AA2=spdiags(centdiag,0,I*I,I*I)+spdiags([updiag;0],1,I*I,I*I)+spdiags([lowdiag;0],-1,I*I,I*I);
    
    
    %if max(abs(sum(A,2)))>10^(-12)
    %    disp('Improper Transition Matrix')
    %    break
    %end
    
    B = (rho + 1/Delta)*speye(I*I) - (AA1 + AA2);
    
    b1 = reshape(pi1 + V1/Delta,[I*I,1]);
    b2 = reshape(pi2 + V2/Delta,[I*I,1]);
    V1 = reshape(B\b1,[I,I]); %SOLVE SYSTEM OF EQUATIONS
    V2 = reshape(B\b2,[I,I]); %SOLVE SYSTEM OF EQUATIONS
    
    V1change = V1 - v1;
    V2change = V2 - v2;
    v1 = V1;
    v2 = V2;
    dist(n) = max(max(max(abs(V1change)), max(abs(V2change))));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
toc;

%surf(v2)
%xlabel("c_2")
%ylabel("c_1")
%zlabel("value_2")







