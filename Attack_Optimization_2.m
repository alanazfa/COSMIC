function [out1,out2,out3,out4] = Attack_Optimization_2(LB_P_load,ps)

global zetta
W = 5;

% load('SF.mat')
if ps.goal ==1
    load('A.mat')
    SF = A;
    SF = SF(ps.load_attack);
    f_int = ps.f_int(1);
    loc = ps.shunt_buses_loc(ps.load_attack);


    cvx_begin
    cvx_solver sedumi

    variables delta_P_L(size(loc,1))

        maximize (  W*(f_int + SF*delta_P_L))
        subject to





         delta_P_L >= -zetta.*LB_P_load(loc);
         delta_P_L <= 0;

    %      F = f_int + SF*delta_P_L;






    cvx_end
elseif ps.goal == 2
    load('A2.mat')
    SF = A2;
    SF = SF(ps.load_attack);
    f_int = 38;
    loc = ps.shunt_buses_loc(ps.load_attack);


    cvx_begin
    cvx_solver sedumi

    variables delta_P_L(size(loc,1))

        maximize (  W*(f_int + SF*delta_P_L))
        subject to





         delta_P_L >= -zetta.*LB_P_load(loc);
         delta_P_L <= 0;

    %      F = f_int + SF*delta_P_L;






    cvx_end
    
else
    load('A.mat')
    load('A2.mat')
    SF = A;
    SF = SF(ps.load_attack);
    SF2 = A2;
    SF2 = SF2(ps.load_attack);
    f_int2 = 38;
    f_int = 103;
    loc = ps.shunt_buses_loc(ps.load_attack);
% 
% 
%     cvx_begin 
% 
%     variables delta_P_L(size(loc,1))
% 
%         maximize (  W*(f_int2 + SF2*delta_P_L)  + W*(f_int + SF*delta_P_L))
%         subject to
% 
% 
% 
% 
% 
%          delta_P_L >= -0.28.*LB_P_load(loc);
%          delta_P_L <= 0;
% 
%     %      F = f_int + SF*delta_P_L;
% 
% 
% 
% 
% 
% 
%     cvx_end
    delta_P_L = -zetta.*LB_P_load(loc);
    cvx_optval=1000;
end
2;
if size(delta_P_L,2)~=51
    new_load_change = zeros(51,1);
    new_load_change(ps.load_attack) = delta_P_L;
    delta_P_L = new_load_change;
    
    percent = zeros(51,1);
    new_percent = delta_P_L(ps.load_attack)./LB_P_load(loc);
    percent(ps.load_attack) = new_percent;
end
F = f_int + SF*delta_P_L(ps.load_attack);


out1 = delta_P_L;
out2 = F;
out3 = cvx_optval;
% out4 = delta_P_L(ps.load_attack)./LB_P_load(loc);
out4 = percent;


end