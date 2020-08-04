function [t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd,P3,Temperature,If,It] = read_outfile(fname,ps,opt)
% usage: [t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd,P3] = read_outfile(fname,ps,opt)

n  = size(ps.bus,1);
ng = size(ps.gen,1);
m  = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices(n,ng,m,n_sh,opt);

data = csvread(fname,1);

t = data(:,1);
X = data(:,1+(1:ix.nx));
Y = data(:,1+ix.nx+(1:ix.ny));

% X vars
delta = X(:,ix.x.delta);
omega = X(:,ix.x.omega)*2*pi*ps.frequency;
Pm    = X(:,ix.x.Pm);
Eap   = X(:,ix.x.Eap);
E1    = X(:,ix.x.E1);
Efd   = X(:,ix.x.Efd);
P3   = X(:,ix.x.P3);
Temperature = X(:,ix.x.temp)+20; % the temerature value in X is reference to 20

% Y vars
Vmag  = Y(:,ix.y.Vmag);
theta = Y(:,ix.y.theta);


V_all = Vmag.*exp(j.*theta);

If_all = ps.Yf*V_all'; % branch status is accounted for in Yf
It_all = ps.Yt*V_all'; % branch status is accounted for in Yt
If = If_all';
It = It_all';
