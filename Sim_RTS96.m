%% simulate 39-bus case
clear all; close all; clc; C = psconstants;
% load('record3.mat')
% load('record2.mat')
load('record6.mat')
% for i = 5:12
number = 0;
number2 = 0;
events_counter = 0;
identifier = NaN(49,2);
total_attack = NaN(49,2);
% for i = 1:49
i = 28;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

% simulation time
t_max = 90;
load rts96_2;
% load RTS_Demand;
% % 
% ps.shunt(:,C.sh.P)= RTS_Demand(:,2);
% select data case to simulate

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
for qq=1:3:size(ps.shunt,1)
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

lembda_list=[0.1:0.1:1.2];

lembda=lembda_list(10);
ps.shunt(:,C.sh.P)     = ps.shunt(:,C.sh.P).*lembda;

limit_gen = find(ps.gen(:,C.ge.Pg)*lembda>=ps.gen(:,C.ge.Pmax));
for g=1:size(limit_gen,1)
    ps.gen(limit_gen(g),C.ge.Pg)    =   ps.gen(limit_gen(g),C.ge.Pmax);
    ps.gov(limit_gen(g),C.go.Pmax)  =   ps.gov(limit_gen(g),C.go.Pmax).*1.1;
end

limit_gen_min = find(ps.gen(:,C.ge.Pg)*lembda<ps.gen(:,C.ge.Pmax));
for gg=1:size(limit_gen_min,1)
    ps.gen(limit_gen_min(gg),C.ge.Pg)    =   ps.gen(limit_gen_min(gg),C.ge.Pg).*lembda;
    ps.gov(limit_gen_min(gg),C.go.Pmax)    = ps.gov(limit_gen_min(gg),C.go.Pmax).*lembda.*1.1;
end

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
opt.sim.oc_limit = 1.75;
opt.sim.uvls_tdelay_ini = 0.2;  % 1 sec delay for uvls relay.
opt.sim.ufls_tdelay_ini = 0.2;  % 1 sec delay for ufls relay.
opt.sim.dist_tdelay_ini = 0.5;  % 1 sec delay for dist relay.
opt.sim.temp_tdelay_ini = 0.1;    % 0 sec delay for temp relay.
opt.sim.oc_tdelay_ini = 0;  % 1 sec delay for oc relay.
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

t_delay([ix.re.oc])= opt.sim.oc_tdelay_ini;

t_prev_check = nan(size(ps.relay,1),1);
dist2threshold = inf(size(ix.re.oc,2)*2,1);
state_a = zeros(size(ix.re.oc,2)*2,1);




%%



% load 'attack_load_2'
% load_attack = unique(randi([1 51],1,15));
load_attack = [2 9 15 21 34 35 37 38 39 41 43 48 49];
% load_attack = [1     4     7    17    23    27    33    39    41    42    47    49    50];
% load_attack = attack_load_2(2).new;
% i=1;
% protected = [ ];
% load_attack = [];
% load_attack= [1:10 20 13 30 31 32 49 51];
% load_attack = [record6(i).load_attack];

% remove_protected = ismember(protected,load_attack);
freq1 = 1;
freq2 = 1;
freq3 = 9;

%% build an event matrix

frac_A = NaN(3,1);
frac_A(1,1)=0.3;
frac_A(2,1)=0.3;
frac_A(3,1)=0.24;
[event,k,step,event_1,event_2,event_3] = event_different_frequancy(C,load_attack,t_max,frac_A,freq1,freq2,freq3);
% EC_time = [11:10:t_max-10];
% for i = 1:size(EC_time,2)
%     EC_event =zeros(1,10);
%     EC_event_less_t  = find(event(:,C.ev.time)<EC_time(1,i));
%     EC_event_more_t = find(event(:,C.ev.time)>=EC_time(1,i));
%     EC_event(1,C.ev.time)      =  EC_time(1,i);
%     EC_event(1,C.ev.type)      =  C.ev.em_control;
%     EC_event(1,C.ev.shunt_loc) =  NaN;
% %     EC_event(1,C.ev.quantity)  =  0;
%     EC_event(1,C.ev.change_by) =  1;
%     event = [event(EC_event_less_t,:);EC_event;event(EC_event_more_t,:)];
% 
% end



% [event,k,step] = event_different_frequancy_2(C,load_attack,t_max,i);
% event_x = event(2:11,:);
% event = [event(1,:)];
2;
% record6(i).event(2:end-1,8)=(record6(i).event(2:end-1,8)./0.9)*0.73;
% record6(i).event(2:end-1,1)= round(record6(i).event(2:end-1,1),4);
% event(2:end-1,1)= round(event(2:end-1,1),4);
% record6(i).event([48:67],:)=[];
% record6(i).event([43 44],:)=[];
% event([19 37  ],:)=[];
% event(find(event(:,6)==6|event(:,6)==40),:)=[];
% event(find(event(:,6)==39),:)=[];

% event([22 53 74 95 116 137 158 179 200 221 242 253],:)=[];

% record6(i).event(61,:)=[];
% record6(i).event(119,:)=[];
% record6(i).event(177,:)=[];

% event = record6(i).event;
% load event_23;
load_attack_percent = sum(ps.shunt(load_attack,2))*0.73/sum(ps.shunt(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

events_counter = events_counter + 1;


%%
% time_step=[20,10,2.5];
%  time_step=[1.25,(1/0.6)];

%%%%%%%%%%%%%%%%%%%
%  EVENT MATRIX   %
%%%%%%%%%%%%%%%%%%%
% time_step=[1,20,10,5,(1/0.3),(1/0.6),2.5,2];
% freq = 2;
% step =   time_step(freq);
% load_frac=[0.1:0.1:0.9];
% frac = 3;
% lembda_list=[0.1:0.1:1.2];
% load_attack = [1:51];
% 
% n_event = (((60/time_step(freq)))*length(load_attack))+2;
% event_number=(60/time_step(freq));
% if mod(60/time_step(freq),2) ~= 0
%    n_event = (((60/time_step(freq))+1)*length(load_attack))+2;
%    event_number=(60/time_step(freq))+1;
% end
% event = zeros(n_event,C.ev.cols);
% % start
% event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
% % trip branches
% k=1;
% time=0;
% 
% for k=1:event_number
%     if mod(k,2)~=0
%         j=1;
%         for i = (k-1)*length(load_attack)+2 : k*length(load_attack)+1
% %       event(i+1,[C.ev.time C.ev.type]) = [10+jj C.ev.trip_branch];
%             event(i,[C.ev.time C.ev.type]) = [10+time C.ev.shed_load];
%             event(i,C.ev.shunt_loc) = load_attack(j);
% %             event(i,C.ev.quantity) = ps_int.shunt(load_attack(j),2)*load_frac(frac);
%             event(i,C.ev.quantity) = load_frac(frac);
%             event(i,C.ev.change_by) = 1;
%             j=j+1;
%         end
%         j=1;
%     else
%         for i = (k-1)*length(load_attack)+2 : k*length(load_attack)+1
% %     event(i+1,[C.ev.time C.ev.type]) = [10+jj C.ev.trip_branch];
%             event(i,[C.ev.time C.ev.type]) = [10+time C.ev.shed_load];
%             event(i,C.ev.shunt_loc) = load_attack(j);
% %             event(i,C.ev.quantity) = -ps_int.shunt(load_attack(j),2)*load_frac(frac);
%             event(i,C.ev.quantity) = -load_frac(frac);
%             event(i,C.ev.change_by) = 1;
%             j=j+1;
%     
%         end
%     end
%     time=time+time_step(freq);
% end
% 
% event(end,[C.ev.time C.ev.type]) = [t_max C.ev.finish];
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%% run the simulation
blockout=0;
k=1;
% event([2 13 164 175 326 337 488 489],:)=[];                
% event([53:613],:)=[];
% [outputs,ps] = simgrid(ps,event,'sim_case39',opt,k,blockout,step,event_1,event_2,event_3,load_attack);
[outputs,ps,blockout] = simgrid_EC(ps,event,'sim_caseRTS_96_EC',opt,k,blockout,step,event_1,event_2,event_3,load_attack);

% if outputs.demand_lost>8500 
%     number = number+1;
%     identifier(i,1) = 1; 
%     total_attack(i,1) = sum(ps.shunt(load_attack,2))*0.8;
%     rec(i).outputs=outputs;
% end
% total_attack(i,2) = outputs.demand_lost;

% [outputs1] = simgrid(ps,record6(i).event,'sim_case39',opt,k,blockout,step);
% if outputs1.demand_lost>8500 
%     number2 = number2+1;
%     identifier(i,2) = 1; 
%     total_attack(i,2) = sum(ps.shunt(record6(i).load_attack,2))*0.8;
% end
% end
%% print the results
fname = outputs.outfilename;
[t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd,P3,Temperature] = read_outfile(fname,ps,opt);
omega_0 = 2*pi*ps.frequency;
omega_pu = omega / omega_0;
%%
figure(1); clf; hold on; 
nl = size(omega_pu,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xtick',[0 600 1200 1800],...
%     'Xlim',[0 50],'Ylim',[0.995 1.008]);
plot(t,omega_pu);
set(gca,'FontSize',16)
ylabel('\omega (pu)','FontSize',22);
xlabel('time (sec.)','FontSize',22);
% PrintStr = sprintf('OmegaPu_P_%s_%s_%s',CaseName, Contingency, Control);
% print('-dpng','-r600',PrintStr)
%%
figure(2); clf; hold on; 
nl = size(theta,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[-0.2 0.5]);
plot(t,theta);
ylabel('\theta','FontSize',18);
xlabel('time (sec.)','FontSize',18);

%%
figure(3); clf; hold on; 
nl = size(Vmag,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,Vmag);
ylabel('|V|','FontSize',18);
xlabel('time (sec.)','FontSize',18);
%%
figure; clf; hold on; 
nl = size(Vmag,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,Vmag(:,1:24));
ylabel('|V|','FontSize',18);
xlabel('time (sec.)','FontSize',18);
figure; clf; hold on; 
nl = size(Vmag,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,Vmag(:,49:73));
ylabel('|V|','FontSize',18);
xlabel('time (sec.)','FontSize',18);

%%
figure(5); clf; hold on; 
nl = size(Pm,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,Pm(:,1:20));
ylabel('Pm','FontSize',18);
xlabel('time (sec.)','FontSize',18);
%%
figure(6); clf; hold on; 
nl = size(delta,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
% plot(t',delta'.*180./pi);
plot(t,delta);
ylabel('Delta','FontSize',18);
xlabel('time (sec.)','FontSize',18);

% figure(7); clf; hold on; 
% nl = size(Eap,2); colorset = varycolor(nl);
% % set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
% plot(t,Eap);
% ylabel('Eap','FontSize',18);
% xlabel('time (sec.)','FontSize',18);
% 
% figure(8); clf; hold on; 
% nl = size(E1,2); colorset = varycolor(nl);
% % set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
% plot(t,E1);
% ylabel('E1','FontSize',18);
% xlabel('time (sec.)','FontSize',18);
% 
% figure(9); clf; hold on; 
% nl = size(Efd,2); colorset = varycolor(nl);
% % set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
% plot(t,Efd);
% ylabel('Efd','FontSize',18);
% xlabel('time (sec.)','FontSize',18);    

figure(10); clf; hold on; 
nl = size(Temperature,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,Temperature);
ylabel('Temperature ( ^{\circ}C)','Interpreter','tex','FontSize',18);
xlabel('time (sec.)','FontSize',18);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Lines analysis         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(11); clf; hold on; 
nl = size(Temperature,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,Temperature(:,[12 24 41 118 119]));
legend({'L:107-203','L:113-215','L:123-217','L:325-121','L:318-223'},'Location','northwest')
ylabel('Temperature Tie Lines ( ^{\circ}C)','Interpreter','tex','FontSize',18);
xlabel('time (sec.)','FontSize',18);

%%
figure(11); clf; hold on; 
nl = size(Temperature,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,Temperature(:,[80:120]));% legend({'L:107-203','L:113-215','L:123-217','L:325-121','L:318-223'},'Location','northwest')
ylabel('Events ','Interpreter','tex','FontSize',24);
xlabel('time (sec.)','FontSize',24);
set(gca,'FontSize',18)

   %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Frequency analysis       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; clf; hold on; 
nl = size(omega_pu,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xtick',[0 600 1200 1800],...
%     'Xlim',[0 50],'Ylim',[0.995 1.008]);
plot(t,omega_pu(:,1:10));
ylabel('\omega (pu)','FontSize',18);
xlabel('time (sec.)','FontSize',18);

%%
figure; clf; hold on; 
nl = size(omega_pu,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xtick',[0 600 1200 1800],...
%     'Xlim',[0 50],'Ylim',[0.995 1.008]);
plot(t,omega_pu(:,11:20));
ylabel('\omega (pu)','FontSize',18);
xlabel('time (sec.)','FontSize',18);

%%
figure; clf; hold on; 
nl = size(omega_pu,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xtick',[0 600 1200 1800],...
%     'Xlim',[0 50],'Ylim',[0.995 1.008]);
plot(t,omega_pu(:,21:30));
ylabel('\omega (pu)','FontSize',18);
xlabel('time (sec.)','FontSize',18);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theta analysis       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11); clf; hold on; 
nl = size(theta,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[-0.2 0.5]);
% plot(t,theta(:,1:10));
% plot(t,theta(:,1),'-x',t,theta(:,2),'--o',t,theta(:,3),':',t,theta(:,4),'-d',t,theta(:,5),'-.')
% plot(t,theta(:,11),'-x',t,theta(:,12),'--o',t,theta(:,13),':',t,theta(:,14),'-d',t,theta(:,15),'-.',t,theta(:,16),'-s',t,theta(:,17),'-o')
plot(t,theta(:,1),'-.o')
hold on
plot(t,theta(:,3),'-.*')
% legend({'bus1','bus2','bus3','bus4','bus5','bus6','bus7','bus8','bus9','bus10'},'Location','northwest','NumColumns',2)
% legend({'bus1','bus2','bus3','bus4','bus5'},'Location','northwest','NumColumns',2)
% legend({'bus11','bus12','bus13','bus14','bus15','bus16','bus17'},'Location','northwest','NumColumns',2)
legend({'bus1','bus3'},'Location','northwest')

ylabel('\theta Area1','FontSize',18);
xlabel('time (sec.)','FontSize',18);
% hold on
% plot(t,theta(:,18:24))
% hold on
% plot(t,theta(:,73),'*')
% hold on
% plot(t,theta(:,39),'*')
% hold on
% plot(t,theta(:,41),'*')
%%
figure(12); clf; hold on; 
nl = size(theta,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[-0.2 0.5]);
plot(t,theta(:,1:24));
ylabel('\theta Area1','FontSize',18);
xlabel('time (sec.)','FontSize',18);
%%
figure(12); clf; hold on; 
nl = size(theta,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[-0.2 0.5]);
plot(t,theta(:,25:34));
ylabel('\theta Area2','FontSize',18);
xlabel('time (sec.)','FontSize',18);
%%
figure(13); clf; hold on; 
nl = size(theta,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[-0.2 0.5]);
% plot(t,theta(:,49:62));
% plot(t,theta(:,73),'--')
plot(t,theta(:,66),'--')
% plot(t,theta(:,25),'-x',t,theta(:,26),'--o',t,theta(:,27),':',t,theta(:,28),'-d',t,theta(:,29),'-.')

hold on
plot(t,theta(:,21),'*')
% hold on
% plot(t,theta(:,47),'*')
% hold on
% plot(t,theta(:,73),'o') 
legend({'bus223','bus318'},'Location','northwest')
ylabel('\theta Area3','FontSize',18);
xlabel('time (sec.)','FontSize',18);

%%
figure(14); clf; hold on; 
nl = size(theta,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[-0.2 0.5]);
% plot(t,theta(:,49:62));
% plot(t,theta(:,73),'--')
plot(t,theta(:,21),'--')
% plot(t,theta(:,25),'-x',t,theta(:,26),'--o',t,theta(:,27),':',t,theta(:,28),'-d',t,theta(:,29),'-.')

hold on
% plot(t,theta(:,21),'*')
% hold on
% plot(t,theta(:,47),'*')
% hold on
plot(t,theta(:,73),'o') 
legend({'bus223','bus318'},'Location','northwest')
ylabel('\theta Area3','FontSize',18);
xlabel('time (sec.)','FontSize',18);
%%
figure; clf; hold on; 
x=[10 10];
y=[0 1];
plot(x,y,'--b','LineWidth',2)
hold on 
x=[10 11.25];
y=[1 1];
plot(x,y,'--b','LineWidth',2)
hold on 
x=[11.25 11.25];
y=[0 1];
plot(x,y,'--b','LineWidth',2)
hold on
x=[11.25 12.5];
y=[0 0];
plot(x,y,'--b','LineWidth',2)
hold on
x=[12.5 12.5];
y=[0 1];
plot(x,y,'--b','LineWidth',2)
hold on 
x=[12.5 13.75];
y=[1 1];
plot(x,y,'--b','LineWidth',2)
hold on
x=[13.75 13.75];
y=[1 0];
plot(x,y,'--b','LineWidth',2)
hold on
x=[13.75 15];
y=[0 0];
plot(x,y,'--b','LineWidth',2)
hold on
x=[15 15];
y=[0 1];
plot(x,y,'--b','LineWidth',2)
hold on 
x=[15 16.25];
y=[1 1];
plot(x,y,'--b','LineWidth',2)
x=[16.25 16.25];
y=[1 0];
plot(x,y,'--b','LineWidth',2)
hold on 
x=[16.25 17.5];
y=[0 0];
plot(x,y,'--b','LineWidth',2)
hold on 
x=[17.5 17.5];
y=[0 1];
plot(x,y,'--b','LineWidth',2)
hold on 
x=[17.5 18.75];
y=[1 1];
plot(x,y,'--b','LineWidth',2)
hold on 
x=[18.75 18.75];
y=[1 0];
plot(x,y,'--b','LineWidth',2)
xlim([0 20])
hold on

x=[10 10];
y=[0 1];
plot(x,y,':r','LineWidth',2)
hold on 
x=[10 11];
y=[1 1];
plot(x,y,':r','LineWidth',2)
hold on 
x=[11 11];
y=[0 1];
plot(x,y,':r','LineWidth',2)
hold on
x=[11 12];
y=[0 0];
plot(x,y,':r','LineWidth',2)
hold on
x=[12 12];
y=[0 1];
plot(x,y,':r','LineWidth',2)
hold on 
x=[12 13];
y=[1 1];
plot(x,y,':r','LineWidth',2)
hold on
x=[13 13];
y=[1 0];
plot(x,y,':r','LineWidth',2)
hold on
x=[13 14];
y=[0 0];
plot(x,y,':r','LineWidth',2)
hold on
x=[14 14];
y=[0 1];
plot(x,y,':r','LineWidth',2)
hold on 
x=[14 15];
y=[1 1];
plot(x,y,':r','LineWidth',2)
x=[15 15];
y=[1 0];
plot(x,y,':r','LineWidth',2)
hold on 
x=[15 16];
y=[0 0];
plot(x,y,':r','LineWidth',2)
hold on 
x=[16 16];
y=[0 1];
plot(x,y,':r','LineWidth',2)
hold on 
x=[16 17];
y=[1 1];
plot(x,y,':r','LineWidth',2)
hold on 
x=[17 17];
y=[1 0];
plot(x,y,':r','LineWidth',2)
xlim([0 20])
ylabel('Events ','Interpreter','tex','FontSize',24);
xlabel('time (sec.)','FontSize',24);
set(gca,'FontSize',18)

%%



figure; clf; hold on; 
x=[10 10];
y=[0 1];
plot(x,y,'-b','LineWidth',2)
hold on 
x=[10 11.25];
y=[1 1];
plot(x,y,'-b','LineWidth',2)
hold on 
x=[11.25 11.25];
y=[0 1];
plot(x,y,'-b','LineWidth',2)
hold on
x=[11.25 12.5];
y=[0 0];
plot(x,y,'-b','LineWidth',2)
hold on
x=[12.5 12.5];
y=[0 1];
plot(x,y,'-b','LineWidth',2)
hold on 
x=[12.5 13.75];
y=[1 1];
plot(x,y,'-b','LineWidth',2)
hold on
x=[13.75 13.75];
y=[1 0];
plot(x,y,'-b','LineWidth',2)
hold on
x=[13.75 15];
y=[0 0];
plot(x,y,'-b','LineWidth',2)
hold on
x=[15 15];
y=[0 1];
plot(x,y,'-b','LineWidth',2)
hold on 
x=[15 16.25];
y=[1 1];
plot(x,y,'-b','LineWidth',2)
x=[16.25 16.25];
y=[1 0];
plot(x,y,'-b','LineWidth',2)
hold on 
x=[16.25 17.5];
y=[0 0];
plot(x,y,'-b','LineWidth',2)
hold on 
x=[17.5 17.5];
y=[0 1];
plot(x,y,'-b','LineWidth',2)
hold on 
x=[17.5 18.75];
y=[1 1];
plot(x,y,'-b','LineWidth',2)
hold on 
x=[18.75 18.75];
y=[1 0];
plot(x,y,'-b','LineWidth',2)
hold on



x=[18.75 20];
y=[0 0];
plot(x,y,'-b','LineWidth',2)
hold on 
x=[20 20];
y=[0 1];
plot(x,y,'-b','LineWidth',2)
hold on 
x=[20 21.25];
y=[1 1];
plot(x,y,'-b','LineWidth',2)
hold on
x=[21.25 21.25];
y=[1 0];
plot(x,y,'-b','LineWidth',2)
hold on
x=[21.25 22.5];
y=[0 0];
plot(x,y,'-b','LineWidth',2)
hold on
x=[22.5 22.5];
y=[0 1];
plot(x,y,'-b','LineWidth',2)
hold on
x=[22.5 23.75];
y=[1 1];
plot(x,y,'-b','LineWidth',2)
hold on
x=[23.75 23.75];
y=[1 0];
plot(x,y,'-b','LineWidth',2)
hold on
x=[23.75 25];
y=[0 0];
plot(x,y,'-b','LineWidth',2)
hold on
x=[25 25];
y=[0 1];
plot(x,y,'-b','LineWidth',2)
hold on
x=[25 26.25];
y=[1 1];
plot(x,y,'-b','LineWidth',2)
hold on
x=[26.25 26.25];
y=[1 0];
plot(x,y,'-b','LineWidth',2)
hold on
x=[26.25 27.5];
y=[0 0];
plot(x,y,'-b','LineWidth',2)
hold on
x=[27.5 27.5];
y=[0 1];
plot(x,y,'-b','LineWidth',2)
hold on
x=[27.5 28.75];
y=[1 1];
plot(x,y,'-b','LineWidth',2)
hold on
x=[28.75 28.75];
y=[1 0];
plot(x,y,'-b','LineWidth',2)
hold on
x=[28.75 30];
y=[0 0];
plot(x,y,'-b','LineWidth',2)
hold on
x=[30 30];
y=[0 1];
plot(x,y,'-b','LineWidth',2)
hold on
x=[30 31.25];
y=[1 1];
plot(x,y,'-b','LineWidth',2)
hold on
x=[31.25 31.25];
y=[1 0];
plot(x,y,'-b','LineWidth',2)
hold on
x=[31.25 32.5];
y=[0 0];
plot(x,y,'-b','LineWidth',2)
hold on
x=[32.5 32.5];
y=[0 1];
plot(x,y,'-b','LineWidth',2)
hold on
x=[32.5 33.75];
y=[1 1];
plot(x,y,'-b','LineWidth',2)
hold on
x=[33.75 33.75];
y=[1 0];
plot(x,y,'-b','LineWidth',2)
hold on

xlim([0 40])
hold on

ylabel('Events ','Interpreter','tex','FontSize',24);
xlabel('time (sec.)','FontSize',24);
set(gca,'FontSize',18)

%%
x=[10 10];
y=[0 1];
plot(x,y,':r','LineWidth',2)
hold on 
x=[10 30];
y=[1 1];
plot(x,y,':r','LineWidth',2)
hold on 
x=[30 30];
y=[0 1];
plot(x,y,':r','LineWidth',2)
hold on
x=[30 50];
y=[0 0];
plot(x,y,':r','LineWidth',2)

ylabel('Events ','Interpreter','tex','FontSize',24);
xlabel('time (sec.)','FontSize',24);
set(gca,'FontSize',18)
