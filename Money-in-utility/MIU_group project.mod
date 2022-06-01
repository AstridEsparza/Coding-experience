% Define Endogeneous Variables


var c, m, n, k, R, rr, pi, z, y, kappa;


% c  for consumption
% m  for real money holding
% n  for labour supply
% k  for production capital
% R  for gross nom. interest rate
% rr for gross real interest rate
% pi for gross inflation rate
% z  for technology shock
% y  for real output
% kappa for monetary policy shock


% Define Exogeneous Variables


varexo varkappa, xi;


% varkappa is the shock of monetary policy
% xi is the shock of technology


% Define Parameters (including Steady-State variables)


parameters beta sigma psi theta_1 theta_2 omega eta rho_z rho_pi delta n_ss c_ss y_ss rr_ss R_ss pi_ss k_ss m_ss;


% beta discount factor
% sigma consumption factor in utility fnc.
% psi money factor in utility fnc.
% theta_1 money factor in utility fnc. 
% theta_2 labour factor in utility fnc.
% omega labour factor in utility fnc.
% eta prod. fnc. share of capital
% rho_z Technology shock autocorrelation - AR(1) Process
% rho_pi monetary policy CB's reaction to inflation 
% delta capital depreciation
% tau monetary injection
% n_ss ss_labour
% c_ss ss_consumption
% y_ss ss_production
% rr_ss ss_real_interest_rate
% R_ss ss_nom_interest_rate
% pi_ss ss_inflation
% k_ss ss_capital
% m_ss ss_money


% Parameter Calibration


eta=0.3;
delta=0.025;
beta=0.99;
sigma=2;
psi=2;
omega=1;
rho_pi=1.5;
rho_z=0.9;
theta_1=12.1873;
theta_2=0.0019;


% steady-state variables


n_ss=1/3;
rr_ss=1/beta;
k_ss=(1/3)*(((1/beta)-1+delta)/eta)^(1/(eta-1));
y_ss=(1/3)*(((1/beta)-1+delta)/eta)^(eta/(eta-1));
c_ss=y_ss-delta*k_ss;
R_ss=(1/beta)^(rho_pi/(rho_pi-1 ));
pi_ss=(1/beta)^(1/(rho_pi-1));
m_ss=((theta_2*((c_ss)^(sigma))*R_ss)/(R_ss-1))^(1/psi);


% Define Linear Model


model(linear);


c*c_ss+k*k_ss=y_ss*y+(1-delta)*k_ss*k(-1); % Resource constraint
y=z+eta*k(-1)+(1-eta)*n; % Production function
-sigma*c=rr-sigma*c(+1); % Euler Equation
rr_ss*rr=(eta*(y_ss/k_ss))*(y(+1)-k); % Real Return on Capital 
sigma*c=(R/(R_ss-1))+psi*m; % Intratemp. Money/Cons. Trade-Off
n=(y-sigma*c)/(1+omega); % Intratemp. Labout/Cons. Trade-Off
R=rr+pi(+1); % Fisher Equation
R=kappa+rho_pi*pi; % Monetary policy rule
%Shocks
z=rho_z*z(-1)+xi;
kappa=varkappa;
end;


shocks;
var xi; stderr 1;
var varkappa; stderr 1;
end;


% Start Simulation
stoch_simul;
