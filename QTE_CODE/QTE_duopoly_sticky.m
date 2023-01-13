% Ivan Rendo Barreiro - QTE - Sticky duopoly (Meza and Rios)
% Modication of Ben Moll's HJB_no_uncertainty_implicit.m code.

tic; 
clear all; clc;


s=0.2; % speed of adjustement
c=1;
rho=0.05; 
a = 10; 

I=4000;
pmin = c+0.01;
pmax = 10;
p = linspace(pmin,pmax,I)';
dp = (pmax-pmin)/(I-1);

maxit=1000;
crit = 10^(-6);

Delta = 1000; % mas grande, mas lento pero mas estable

dV1f = zeros(I,1);
dV2f = zeros(I,1);
dV1b = zeros(I,1);
dV2b = zeros(I,1);
q1 = p;
q2 = p;

%INITIAL GUESS
v0 = ((p-c-(a-p)./2).*((a-p)./2))./rho; % net present value of profits in the SS, assuming q1=q2
v1 = v0;
v2 = v0;

for n=1:maxit
    V1 = v1;
    V2 = v2;
    % forward difference
    dV1f(1:I-1) = (V1(2:I)-V1(1:I-1))/dp;
    dV2f(1:I-1) = (V2(2:I)-V2(1:I-1))/dp;
    dV1f(I) = (p(I)-c-0)/s ; % using FOC and q>=0
    dV2f(I) = (p(I)-c-0)/s; % using FOC and q>=0
    % backward difference
    dV1b(2:I) = (V1(2:I)-V1(1:I-1))/dp;
    dV2b(2:I) = (V2(2:I)-V2(1:I-1))/dp;
    dV1b(1) = (p(1)-c-0)/s; % using FOC and q>=0
    dV2b(1) = (p(1)-c-0)/s; % using FOC and q>=0
    
    % quantities and price change with forward difference
    q1f = max(p - c - s.*dV1f,0);
    q2f = max(p - c - s.*dV2f,0);
    pf = s.*(a - q1f - q2f - p); % by definition of p'
    % quantities and price change backward difference
    q1b = max(p - c - s.*dV1b,0);
    q2b = max(p - c - s.*dV2b,0);
    pb = s.*(a - q1b - q2b - p); % by definition of p'
    
    %consumption and derivative of value function at steady state
    con = 0;
    dV10 = zeros(I,1) + con; %(p-c-a+q2 + p)./s; This is a trick...
    dV20 = zeros(I,1) + con; %(p-c-a+q1 + p)./s; 
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = pf > 0; %positive drift --> forward difference
    Ib = pb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); 
    
    dV1_Upwind = dV1f.*If + dV1b.*Ib + dV10.*I0;
    dV2_Upwind = dV2f.*If + dV2b.*Ib + dV20.*I0;
    
    q1 = max(p - c - s.*dV1_Upwind,0);
    q2 = max(p - c - s.*dV2_Upwind,0);
    pi1 = (q1.*(p-c));
    pi2 = (q2.*(p-c));

    %CONSTRUCT MATRIX
    X = -min(pb,0)/dp;
    Y = -max(pf,0)/dp + min(pb,0)/dp;
    Z =  max(pf,0)/dp;
    
    A = spdiags(Y,0,I,I)+spdiags(X(2:I),-1,I,I)+spdiags([0;Z(1:I-1)],1,I,I);
    
    if max(abs(sum(A,2)))>10^(-12)
        disp('Improper Transition Matrix')
        break
    end
    
    B = (rho + 1/Delta)*speye(I) - A;
    
    b1 = pi1 + V1/Delta;
    b2 = pi2 + V2/Delta;
    V1 = B\b1; %SOLVE SYSTEM OF EQUATIONS
    V2 = B\b2; %SOLVE SYSTEM OF EQUATIONS
    V1change = V1 - v1;
    V2change = V2 - v2;
    v1 = V1;
    v2 = V2;
    dist(n) = max(max(abs(V1change)), max(abs(V2change)));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
%plot(p,q1)
toc;

pdot = s.*(a-q1-q2-p);
T = 1000;
price = zeros(1,T);
price(1) = 6.2;
delta=0.0008;
for i=2:T
    [a,b] = min(abs(p-price(i-1)));
    price(i) = price(i-1) + pdot(b)*delta;
end

plot(1:T,price)
hold on;

T = 1000;
price = zeros(1,T);
price(1) = 6;
delta=0.0008;
for i=2:T
    [a,b] = min(abs(p-price(i-1)));
    price(i) = price(i-1) + pdot(b)*delta;
end
plot(1:T,price)

T = 1000;
price = zeros(1,T);
price(1) = 6.5;
delta=0.0008;
for i=2:T
    [a,b] = min(abs(p-price(i-1)));
    price(i) = price(i-1) + pdot(b)*delta;
end
plot(1:T,price)
hold off

plot(p,q1)






