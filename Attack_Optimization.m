function [out1,out2,out3] = Attack_Optimization(LB_P_load,ps)

C = -ones(1,11);
W = 1;
M = 2000;
load('SF.mat')
SF = SF(ps.load_attack);
f_int = ps.f_int;
loc = ps.shunt_buses_loc(ps.load_attack);
addpath 'C:\Users\alanazfa\OneDrive\CVX\HW\cvx\9.1\toolbox\R2015a'
addpath 'C:\Users\alanazfa\OneDrive\CVX\HW\cvx\9.1\tools\platform\win64x86\bin'
addpath 'C:\Users\alanazfa\OneDrive\CVX\HW\cvx\mosek\w64'

cvx_begin 
cvx_solver mosek
variables delta_P_L(11) F(1)
variable B integer

    minimize (   C*delta_P_L + W*F)
    subject to
            

     
    

     delta_P_L >= -0.3.*LB_P_load(loc);
     delta_P_L <= 0;
     
     (f_int + SF*delta_P_L)+ M*B - F >= 0;
     (f_int + SF*delta_P_L)+ M*B + F <= M;
     (f_int + SF*delta_P_L) <= F;
     -(f_int + SF*delta_P_L) <= F;
     B <= 1;
     
  

     
  
cvx_end
2;



out1 = delta_P_L;
out2 = F;
out3 = cvx_optval;




end