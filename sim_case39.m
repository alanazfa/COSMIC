%% simulate 39-bus case
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

% simulation time
t_max = 130;

% select data case to simulate
ps = updateps(case39_ps);
ps = updateps(ps);
% ps = replicate_case(ps,2);          
ps = unify_generators(ps); 
ps.branch(:,C.br.tap)       = 1;
ps.shunt(:,C.sh.factor)     = 1;
ps.shunt(:,C.sh.status)     = 1;
% ps.shunt(:,C.sh.frac_S)     = 0;
% ps.shunt(:,C.sh.frac_E)     = 0;
% ps.shunt(:,C.sh.frac_Z)     = 1;
% ps.shunt(:,C.sh.gamma)      = 0.08;

[ZIP_Coff_P_Q_Final] = ZIP_Coeff();
for qq=3:3:size(ps.shunt,1)
    % real power for residential
    ps.shunt(qq,C.sh.frac_S)     = ZIP_Coff_P_Q_Final.Pp_2(1);
    ps.shunt(qq,C.sh.frac_Z)     = ZIP_Coff_P_Q_Final.Zp_2(1);
    % reactive power for residential
    ps.shunt(qq,C.sh.frac_Q_Z)     = ZIP_Coff_P_Q_Final.Zq_2(1);
    ps.shunt(qq,C.sh.frac_Q_I)     = ZIP_Coff_P_Q_Final.Iq_2(1);
    ps.shunt(qq,C.sh.frac_Q_S)     = ZIP_Coff_P_Q_Final.Pq_2(1);
    
    if qq+1<=size(ps.shunt,1)
    % real power for Commercial
        ps.shunt(qq+1,C.sh.frac_S)     = ZIP_Coff_P_Q_Final.Pp_2(2);
        ps.shunt(qq+1,C.sh.frac_Z)     = ZIP_Coff_P_Q_Final.Zp_2(2);
    % reactive power for Commercial
        ps.shunt(qq+1,C.sh.frac_Q_Z)     = ZIP_Coff_P_Q_Final.Zq_2(2);
        ps.shunt(qq+1,C.sh.frac_Q_I)     = ZIP_Coff_P_Q_Final.Iq_2(2);
        ps.shunt(qq+1,C.sh.frac_Q_S)     = ZIP_Coff_P_Q_Final.Pq_2(2);
    end
    % real power for Industrial
     if qq+2<=size(ps.shunt,1)
        ps.shunt(qq+2,C.sh.frac_S)     = ZIP_Coff_P_Q_Final.Pp_2(3);
        ps.shunt(qq+2,C.sh.frac_Z)     = ZIP_Coff_P_Q_Final.Zp_2(3);
    % reactive power for Industrial
        ps.shunt(qq+2,C.sh.frac_Q_Z)     = ZIP_Coff_P_Q_Final.Zq_2(3);
        ps.shunt(qq+2,C.sh.frac_Q_I)     = ZIP_Coff_P_Q_Final.Iq_2(3);
        ps.shunt(qq+2,C.sh.frac_Q_S)     = ZIP_Coff_P_Q_Final.Pq_2(3);
     end
   
end
% ps = replicate_case(ps,2);          
ps = unify_generators(ps); 
ps.branch(:,C.br.tap)       = 1;
ps.shunt(:,C.sh.factor)     = 1;
ps.shunt(:,C.sh.status)     = 1;
ps.shunt(:,C.sh.frac_S)     = 1;
ps.shunt(:,C.sh.frac_E)     = 0;
ps.shunt(:,C.sh.frac_Z)     = 0;
ps.shunt(:,C.sh.gamma)      = 0.08;

% to differentiate the line MVA ratings
rateB_rateA                     = ps.branch(:,C.br.rateB)./ps.branch(:,C.br.rateA);
rateC_rateA                     = ps.branch(:,C.br.rateC)./ps.branch(:,C.br.rateA);
ps.branch(rateB_rateA==1,C.br.rateB)    = 1.1 * ps.branch(rateB_rateA==1,C.br.rateA);
ps.branch(rateC_rateA==1,C.br.rateC)    = 1.5 * ps.branch(rateC_rateA==1,C.br.rateA);

% set some options
opt = psoptions;
opt.sim.integration_scheme = 1;
opt.sim.dt_default = 1/10;
opt.nr.use_fsolve = true;
% opt.pf.linesearch = 'cubic_spline';
opt.verbose = true;
opt.sim.gen_control = 1;        % 0 = generator without exciter and governor, 1 = generator with exciter and governor
opt.sim.angle_ref = 0;          % 0 = delta_sys, 1 = center of inertia---delta_coi
                                % Center of inertia doesn't work when having islanding
opt.sim.COI_weight = 0;         % 1 = machine inertia, 0 = machine MVA base(Powerworld)
opt.sim.uvls_tdelay_ini = 0.5;  % 1 sec delay for uvls relay.
opt.sim.ufls_tdelay_ini = 0.5;  % 1 sec delay for ufls relay.
opt.sim.dist_tdelay_ini = 0.5;  % 1 sec delay for dist relay.
opt.sim.temp_tdelay_ini = 0;    % 0 sec delay for temp relay.
% Don't forget to change this value (opt.sim.time_delay_ini) in solve_dae.m

% initialize the case
ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
ps = update_load_freq_source(ps);
% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
% initialize relays
ps.relay                    = get_relays(ps,'all',opt);

% initialize global variables
global t_delay t_prev_check dist2threshold state_a 
n    = size(ps.bus,1);
ng   = size(ps.mac,1);
m    = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices(n,ng,m,n_sh,opt);
t_delay = inf(size(ps.relay,1),1);
t_delay([ix.re.uvls])= opt.sim.uvls_tdelay_ini;
t_delay([ix.re.ufls])= opt.sim.ufls_tdelay_ini;
t_delay([ix.re.dist])= opt.sim.dist_tdelay_ini;
t_delay([ix.re.temp])= opt.sim.temp_tdelay_ini;
t_prev_check = nan(size(ps.relay,1),1);
dist2threshold = inf(size(ix.re.oc,2)*2,1);
state_a = zeros(size(ix.re.oc,2)*2,1);

%% build an event matrix
load_attack = 1:19
% load_attack_area_1 = [2 3 4 5 18 19];
load_attack_area_1 = [2 ];
load_attack_area_2 = [6 7 9 10 11 12];
% load_attack_area_3 = [1 8 13 14 15 16 17];
load_attack_area_3 = [1 8 15 ];
% freq1 = rec_events(itr).freq(1);
% freq2 = rec_events(itr).freq(2);
% freq3 = rec_events(itr).freq(3);
freq1 = 1;
freq2 = 9;
freq3 = 1;
frac_A(1,1)=0.3;
frac_A(2,1)=0.3;
frac_A(3,1)=0.3;
[event,k,step,event_1,event_2,event_3] = event_different_frequancy_39(C,load_attack,t_max,frac_A,freq1,freq2,freq3,load_attack_area_1,load_attack_area_2,load_attack_area_3);
% area2_loads = find(event(:,6)>17 & event(:,6)<35);
% event(area2_loads,:)=[];

% event = zeros(4,C.ev.cols);
% % start
% event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
% % trip a branch
% event(2,[C.ev.time C.ev.type]) = [3 C.ev.trip_branch];
% event(2,C.ev.branch_loc) = 35;
% % trip a branch
% event(3,[C.ev.time C.ev.type]) = [3 C.ev.trip_branch];
% event(3,C.ev.branch_loc) = 23;
% % set the end time
% event(4,[C.ev.time C.ev.type]) = [t_max C.ev.finish];
blockout =0;
%% run the simulation
[outputs,ps] = simgrid_attack_toplogy(ps,event,'sim_caseRTS_96_EC',opt,k,blockout,step,event_1,event_2,event_3,load_attack);

%% print the results
fname = outputs.outfilename;
[t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd,P3,Temperature] = read_outfile(fname,ps,opt);
omega_0 = 2*pi*ps.frequency;
omega_pu = omega / omega_0;

figure(1); clf; hold on; 
nl = size(omega_pu,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xtick',[0 600 1200 1800],...
%     'Xlim',[0 50],'Ylim',[0.995 1.008]);
plot(t,omega_pu);
ylabel('\omega (pu)','FontSize',18);
xlabel('time (sec.)','FontSize',18);
% PrintStr = sprintf('OmegaPu_P_%s_%s_%s',CaseName, Contingency, Control);
% print('-dpng','-r600',PrintStr)

figure(2); clf; hold on; 
nl = size(theta,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[-0.2 0.5]);
plot(t,theta);
ylabel('\theta','FontSize',18);
xlabel('time (sec.)','FontSize',18);


figure(3); clf; hold on; 
nl = size(Vmag,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,Vmag);
ylabel('|V|','FontSize',18);
xlabel('time (sec.)','FontSize',18);


figure(5); clf; hold on; 
nl = size(Pm,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,Pm);
ylabel('Pm','FontSize',18);
xlabel('time (sec.)','FontSize',18);

figure(6); clf; hold on; 
nl = size(delta,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
% plot(t',delta'.*180./pi);
plot(t,delta);
ylabel('Delta','FontSize',18);
xlabel('time (sec.)','FontSize',18);

figure(7); clf; hold on; 
nl = size(Eap,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,Eap);
ylabel('Eap','FontSize',18);
xlabel('time (sec.)','FontSize',18);

figure(8); clf; hold on; 
nl = size(E1,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,E1);
ylabel('E1','FontSize',18);
xlabel('time (sec.)','FontSize',18);

figure(9); clf; hold on; 
nl = size(Efd,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,Efd);
ylabel('Efd','FontSize',18);
xlabel('time (sec.)','FontSize',18);    

figure(10); clf; hold on; 
nl = size(Temperature,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,Temperature);
ylabel('Temperature ( ^{\circ}C)','Interpreter','tex','FontSize',18);
xlabel('time (sec.)','FontSize',18);


