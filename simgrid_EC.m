function [outputs,ps,blockout] = simgrid_EC(ps,event,outfilename,opt,k,blockout,step,event_1,event_2,event_3,load_attack)
% usage: [outputs,ps] = simgrid(ps,event,outfilename,opt)
% perform a cascading failure simulation
% inputs:
%  ps           - powersystem structure. see psconstants.
%  event        - a matrix defining the exogenous events to simulate. see psconstants.
%  outfilename  - the name for the output file
%  opt          - an option structure. see psoptions.
% outputs:
%  outputs.         - a structure of output data, including the following:
%    .success       - a bool flag indicating whether we were able to simulate to the end.
%    .t_simulated   - vector of time periods in the final results.
%    .outfilename   - the name of the file where the results are stored.
%    .demand_lost   - the amount of demand that was shed during the simulation.
%    .computer_time - the amount of wall-clock time required to perform this simulation
%    .endo_events   - matrix that logs endogenous events during the simulation.
%  ps               - the powersystem structure at the end of the simulation.

%% process the inputs
if nargin<2 || isempty(event)
    if ~isempty(ps.event) 
        event = ps.event; 
    else
        error('simgrid:err','undefined event matrix.');
    end
    if ~isfield('branch_i')
        error('simgrid:err','need ps.branch_i to run simgrid');
    end
end
if nargin<3 || isempty(outfilename)
    outfilename = 'outfile.csv';
end
if nargin<4 || isempty(opt)
    opt = psoptions;
end

%% clean up the ps structure
ps = updateps(ps);

%% edit the output file name
start_time = datestr(now,30);
outfilename = sprintf('%s_%s.csv',strtok(outfilename,'.'),start_time);
tracefilename = sprintf('trace_%s_%s.mat',strtok(outfilename,'.'),start_time);

%% prepare the outputs
outputs.success         = false;
outputs.t_simulated     = [];
outputs.outfilename     = outfilename;
outputs.demand_lost     = 0;
outputs.computer_time   = [];   ct = clock;
outputs.start_time      = start_time;     
ps.event_record         = [];
event_record            = [];
outputs.linear_solves   = [];

%% other prep work
C           = psconstants;
verbose     = opt.verbose;
discrete    = false;

%% check the event matrix
event = sortrows(event);
if event(1,C.ev.type)~=C.ev.start
    error('simgrid:err','first event should be a start event.');
end
t_0   = event(1,1);
event = event(2:end,:); % remove the start event
if event(end,C.ev.type)~=C.ev.finish
    error('simgrid:err','last event should be a finish event.');
end
t_end = event(end,1);
t = t_0;

%% print something
out = fopen(outfilename,'w');
if isempty(out)
    error('simgrid:err','could not open outfile: %s',outfilename);
end
fprintf(out,'starting simulation at t = %g\n',t);
fclose(out);

if verbose
    fprintf('starting simulation at t = %g\n',t);
    fprintf('writing results to %s\n',outfilename);
end

%% get Ybus, x, y
% build the Ybus
if ~isfield(ps,'Ybus') || ~isfield(ps,'Yf') || ~isfield(ps,'Yt')
    [ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
end
% build x and y
[x,y] = get_xy(ps,opt);

%% step through the simulation
event_counter=0;
recoreded = 0;
attack_index1 = 0;
attack_index2 = 0;
attack_index3 = 0;
attack_index = 0;
load_location=[];
event_no = 1;
load_EV1 = [];
load_EV2 = [];
load_EV3 = [];
Load_buses = ps.shunt(:,[1 2]);
ps.Load_buses = Load_buses;
t_indicator = 0;
ps.load_attack = load_attack;
shunt_buses_loc = NaN(size(ps.shunt(:,1),1),1);
for loc_shunts = 1:size(ps.shunt(:,1),1)
    shunt_buses_loc(loc_shunts,1) = find(ps.shunt(loc_shunts,1)== ps.bus(:,1));
end
ps.shunt_buses_loc = shunt_buses_loc;    
while t < t_end
    % process exogenous discrete events at this time period
    while event(event_no,1) <= t
        [ps,new_discrete] = process_event(ps,event(event_no,:),opt);
        EV = event(event_no,:);
        if EV(C.ev.type)~=22 & EV(C.ev.type)~=23
            sheded_amount = ps.shunt(EV(C.ev.shunt_loc),2)*EV(C.ev.quantity);
            Load_buses(EV(C.ev.shunt_loc),2) = Load_buses(EV(C.ev.shunt_loc),2) - sheded_amount;
            ps.Load_buses(EV(C.ev.shunt_loc),2) = Load_buses(EV(C.ev.shunt_loc),2);
            event_no = event_no+1;
            discrete = discrete | new_discrete;     % for simultaneous events
        else
            event_no = event_no+1;
        end
    end
     
    if discrete
        t = t + opt.sim.t_eps;
    end
    % find the next time interval
    t_next = event(event_no,1);
    % try to simulate between t and t_next
    if opt.verbose
        fprintf('\n Simulating from t=%g s to t=%g s\n',t,t_next);
    end
    [ps,t_out,X,Y,blockout] = simgrid_interval(ps,t,t_next,x,y,opt,blockout);
    
    if opt.verbose  & ~isempty(t_out)
        fprintf(' Completed simulation up until t=%g\n',t_out(end));
    end
    
    % log to files
    if opt.sim.writelog
        if opt.verbose 
            fprintf(' Writing simulation results to %s\n',outfilename);
        end
        write_state(outfilename,t_out,X,Y);
    end
    
    % update time and solutions for next interval
    t = t_next;
    if ~isempty(t_out)
        x = X(:,end);
        y = Y(:,end);
    end

    % record events
    event_record            = [event_record; ps.event_record]; %#ok<AGROW>
    % finding shedded loads by the relays and that load not attacked
    % recording the load to modify it before we solve EC optimization
    % problem
    if ~isempty(ps.event_record)
        if sum(ps.event_record(:,2) == 13) > 0 | sum(ps.event_record(:,2) == 14)>0 | sum(ps.event_record(:,2)==15)>0 ;
            Loads_effected = find(ps.event_record(:,2) ~= 16); 
            for i = 1:size(Loads_effected,1)
                EV_R = ps.event_record(Loads_effected(i),:);
                if sum(load_attack== EV_R(C.ev.shunt_loc)) & (EV_R(C.ev.type)~= 11 | EV_R(C.ev.type)~= 12)
                    sheded_amount = ps.shunt(EV_R(C.ev.shunt_loc),2)*EV_R(C.ev.quantity);
                    Load_buses(EV_R(C.ev.shunt_loc),2) = Load_buses(EV_R(C.ev.shunt_loc),2) - sheded_amount;
                end
            end
        end
    end
        
    ps.event_record         = [];
    if size(find(isnan(x)),1)> 30
        blockout = 2;
        break;
    end
    
    if blockout == 1
        break;
    end
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%   Falah added to adjust event matrix   %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     [event_1,attack_index1,load_EV1] = event1(attack_index1,event_record,event_1,t,load_EV1);
%     [event_2,attack_index2,load_EV2] = event2(attack_index2,event_record,event_2,t,load_EV2);
%     [event_3,attack_index3,load_EV3] = event3(attack_index3,event_record,event_3,t,load_EV3);
%     
% 
% 
%     if size(event_record,1)~=0
%         event_1_load =find(ismember(event(:,6), event_1(:,6) ));
%         event(event_1_load,:) = event_1;
% 
%         event_2_load =find(ismember(event(:,6), event_2(:,6) ));
%         event(event_2_load,:) = event_2;
% 
%         event_3_load =find(ismember(event(:,6), event_3(:,6) ));
%         event(event_3_load,:) = event_3;
%     end

   
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%   Shedding   %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    [P_g_EC,P_L_EC,theta_EC,P_line_EC,cvx_optval_EC,load_locations] = DC_power_flow_RTS_96(load_attack,ps,Load_buses);

%     if ~isempty(event_record) & cvx_optval_EC ~= inf 
%         [EC_load_shedding_loc,event,t,event_1,event_2,event_3,t_indicator] = EC_event(ps,P_L_EC,Load_buses,load_locations,load_attack,event,t,event_1,event_2,event_3,t_indicator,cvx_optval_EC);
%     end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%    event_perecentage_DC = NaN(size(load_attack,2),1);
%    
%    event_perecentage_DC =  ( ps.shunt(load_attack',2) - round(P_L(targeted_load_loc))) ./ps.shunt(load_attack',2);
%    
%    shed_bus = find(event_perecentage_DC~=0);
%    clear_bus = find(event_perecentage_DC == 0);
%    for i = 


end

%% update outputs and clean up
if t_out(end) < (t_end - 1)
    outputs.success         = false;
    outputs.demand_lost     = sum(ps.shunt(:,C.sh.P));
else
    outputs.success         = true;
    outputs.demand_lost     = sum(ps.shunt(:,C.sh.P)) - ps.shunt(:,C.sh.P)'*ps.shunt(:,C.sh.factor);
end
outputs.t_simulated     = [t_0 t_out(end)];
if opt.sim.writelog
    outputs.outfilename     = outfilename;
else 
    outputs.outfilename     = [];
end
outputs.event_record    = event_record;
outputs.computer_time   = etime(clock,ct);

save(tracefilename,'x','y');
if opt.verbose
    fprintf('Completed simulation from %d sec. to %d sec. \n',t_0,t_next);
end