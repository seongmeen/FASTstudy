clear all;
clc;
%% PSD behavior analysis tool
%FAST parameters
%absorber
radius_absorber=0.003;                               r_ab =radius_absorber;                         %[m]
height_absorber=1;                                   h_ab =height_absorber;
thickness_cladding_absorber=0.0001;                  t_ab =thickness_cladding_absorber;
density_absorber=1.2965*(0.85)*(100^3/1000);         den_ab=density_absorber;                     %[kg/m3]
%void and SiC
radius_void=0.003;                                   r_v  =radius_void;
height_void=0.7;                                     h_v  =height_void;
thickness_void=0.0001;                               t_v  =thickness_void;
density_void=0.58721;                                den_v=density_void;                         %[kg/m3] at 545 'C, 0.1Mpa
density_SiC =3.18*(100^3/1000);                      den_SiC=density_SiC;                    
%height and diameter of FASa
H=h_ab+h_v+4*t_ab;                     
%% test
D=2*(r_ab+t_ab);                R=D/2;


  


%% Subchannel parameter


%Pin parameters
daimeter_pin=0.009675;                               D_p   =daimeter_pin;
radius_pin=0.009675/2;                               R_p   =radius_pin;
height_pin=2.9;                                      H_p   =height_pin;    % notion 긴급 코드 수정 확인
height_actual=1.6;                                   H_a   =height_actual; % active core length BOL 일때는 1 m , EOL 일때는 1.6 m 임을 고려
%Volume and Area of FAST
area_side  =pi*(2*(r_ab+t_ab))*(h_ab+h_v+4*t_ab);    A_side =area_side;                          %side area for FAST [m2]
area_front =pi*(r_ab+t_ab)^2;                        A_front=area_front;                         %upper and lower area
area_flow=pi*(R_p^2-(r_ab+t_ab)^2);                  A_flow =area_flow;                          %flow area
volume_FAST=A_front*(h_ab+h_v+4*t_ab);               V =volume_FAST;                             %volume of FAST     [m3]
H_actual=1.6;                                                              %바닥면에서부터 시작한 높이. 

%density of FAST
density_FAST=( (pi*r_ab^2*h_ab)*den_ab+(pi*r_v^2*h_v)*den_v +...
    (V-pi*(r_ab^2*h_ab+r_v^2*h_v))*den_SiC )/ V;  den_FAST=density_FAST;

%hydraulics parameters
hydraulics_diameter=2*R_p-2*(r_ab+t_ab);        d_h=hydraulics_diameter;   %hydraulics diameter De for annular flow
e=8.7E-10;           g=-9.80665;                                           %[m]
gap =R_p-(r_ab+t_ab);                                                      %gap btw FAST and Pin
%% mesh or geometry generation
%meshing generation for FAST
node_MATRA=150;                                                            %MATRA-LMR number of output node
num_h=node_MATRA;                                         dh=H/(num_h-1);
zi_F=0:dh:H;


zi_transient_F=zi_F+(H_p-H);                                               %actual position of FAST per height    
ab_transient_F=H_p-h_ab;                                                   %Absorber height from the bottom
%meshing generation for coolant region
num_c=150;                                                                 %z-direction node number of coolant
num_gap=20;                                                                %r-direction node number of coolant in gap

mesh_cz=H_p/(num_c-1);                             dz=mesh_cz;             %z-direction mesh size of coolant
mesh_cr=(R_p-D/2)/(num_gap-1);                     dr=mesh_cr;             %r-direction mesh size of coolant
ri=D/2:dr:D_p/2;                                                           %r-direction mesh numbering for FAST
dA=A_side/(H/dh);                                                          %differensial side area ( it should be devided by coolant mesh size)


%% Transient temperature calculation for code valdidation assuming no heat flux
%properties of liquid sodium
%k=92.951-5.8087e-2*(T-273.15)+11.7274e-6*(T-273.15)^2                      [W/m.K]
%p=950.0483-0.2298*(T-273.15)-14.6045e-6*(T-273.15)^2+5.6377e-9*(T-273.15)^3[Kg/m3]
%cp=1436.715-0.5806*(T-273.15)+4.6273e-4*(T-273.15)^2                       [J/Kg.K]
%mu=log10(mu)=-2.4892+220.65/T-0.4295*log10(T)                              [Pa s]
%code coupling 때 다시 하는 것으로 한다. MATRA-LMR 위치별 propertives 를 여기로 가져오는 것.
%initial condition
a=0;   v=0;    s=0;                                                        %initial velocity and moving distance of FAST & acceleration
F_drag=0;       F_pressure=0;                                                           %initial Drag force is zero and pressure is 0
dt=0.1;      t=0;                 
%time step and initial time
Pe_dz=0;                                                                   %initial pressure gradient is zero. hydrostatic pressure gradient has been already considered in buoyancy
Q_error=0;     
tol1=1e-5;                                                                  %tolerance        

T=465;                                                                     %constant coolant temperature
T_c_steady=(T+273)*ones(num_c,1);                                          %sodium temperature, 545C outlet , HT9=650 C limitation <여기서 MATRA OUTPUT axial step을 아주 짧게 조절할 수 없다>
k_c=zeros(node_MATRA,1); den_c=zeros(node_MATRA,1);        %size setting for calculation speed
cp_c=zeros(node_MATRA,1); mu_c=zeros(node_MATRA,1);

for i=1:node_MATRA
%liquid sodium properties
    k_c(i,1)=92.951-5.8087e-2*(T_c_steady(i,1)-273.15)+11.7274e-6*(T_c_steady(i,1)-273.15)^2;                                       %k matrix
    den_c(i,1)=950.0483-0.2298*(T_c_steady(i,1)-273.15)-14.6045e-6*(T_c_steady(i,1)-273.15)^2+5.6377e-9*(T_c_steady(i,1)-273.15)^3; %density matrix
    cp_c(i,1)=1436.715-0.5806*(T_c_steady(i,1)-273.15)+4.6273e-4*(T_c_steady(i,1)-273.15)^2;                                        %capacity matrix
    mu_c(i,1)=10^(-2.4892+220.65/T_c_steady(i,1)-0.4295*log10(T_c_steady(i,1)));                                                    %mu matrix   
end



%% Reactivity 
% kinetics parameters
%PKE.beta =[1.66E-04,1.07E-03,1.07E-03,3.26E-03,1.18E-03,3.76E-04];  % BOL
%PKE.ramda=[1.25E-02,3.15E-02,1.11E-01,3.21E-01,1.35E+00,8.87E+00];  % BOL
PKE.beta =[9.42E-05,7.45E-04,6.77E-04,1.96E-03,7.73E-04,2.07E-04]; % EOL
PKE.ramda=[1.25E-02,3.04E-02,1.09E-01,3.21E-01,1.34E+00,9.94E+00]; % EOL
PKE.beta_eff=sum(PKE.beta);
%EOL.beta_eff=sum(EOL.beta);

b=PKE.beta_eff;             %[effective delayed neutron fraction] # 하탄토 논문 참고
dollor=b;
cent=0.01*dollor;
GT=0.768*10^-6;           %prompt neutron generation time[sec]

neg_f=-0.072*cent;         %negative feedback coefficient of the fuel
neg_c=0.201*cent;         %negative feedback coefficient of the coolant
neg_a=-0.06*cent;         %vesel axial expansion
neg_r=-0.153*cent;         %vesel radial expansion 
neg_CEDL=-0.063*cent;      %control rod expansion and contraction
RW=(433.993-756.446)*cent; %Boron-10 20% depletion
RW=RW/3/84;                                                                % 3개씩 84개의 assembly 에 들어있으므로.
Absorber_Reactivity_Ratio=RW*[1,1,1] ;                                     % 맨 왼쪽이 아래쪽임
absorber_section_number=length(Absorber_Reactivity_Ratio);                  % Section number that deviding absorber height

% Initial array
Neutron=zeros(5000,1);
Precursors=[zeros(5000,1),zeros(5000,1),zeros(5000,1),zeros(5000,1),zeros(5000,1),zeros(5000,1)];
% Inital condition
Neutron(1)=1;
for i=1:6
    Precursors(1,i)=PKE.beta(i)/GT/PKE.ramda(i)*Neutron(1);
end


inserted_reactivity_FAST = 0;        %initial inserted reactivity
feedback_reactivity      = 0;        %initial feedback reactivity
rod_withdrawal_reactivity= 0.0048777;%3.5917e-04; % [reactivity]
rod_withdrawal_time=15; % 15second rod withdrawal time

total_reactivity=(inserted_reactivity_FAST)+(feedback_reactivity);


%% stability
stability=-2*mu_c(1,1)/den_c(1,1)*dt/dr^2+1;
%coolant area matrix processing
%% velocity filed.
v_f=zeros(num_gap*num_c,1);                                                %initail velocity field(v=0)  &%size setting for calculation speed
v_f_p=v_f;                                                                 %n step velocity field



   
    
%% FAST region temperature transfer time (lumped capacity method)
%공통 변수
%% interpolation for conductivity
temperature = [ 293 373:100:1173]; %[K] 
weight =5:5:10;   %[%] 

uzr=[15.8000   17.8000   19.8000   21.8000   24.8000   27.4000   30.4000   33.4000   37.0000   40.0000];
S=griddedInterpolant(temperature,uzr);
     
%% Fuel rod side

f_o   =0.0042975;                                                          %fuel outter radius
q=12740;                                                                   %heat generation 392.
con_q=q/(pi*f_o^2*(2*H_a/pi));                                                %0.63662는 cos를 높이로 적분한 값에 비롯됨. 77page 연구노트
D_w=0.0013545;                                                             %wire rap diameter
PD=1.14;                                                                   %P/D ratio

D_fh=4*(3^0.5/4*(D_p*PD)^2-pi*D_p^2/8-pi*D_w^2/8)/(pi/2*(D_p+D_w));        %열 수력학적 직경(m)  %0.004176 - matra
A_c_f=(3^0.5/4*(D_p*PD)^2-pi*D_p^2/8-pi*D_w^2/8);                          %물이 흐르는 면적
%D_fh=0.004176;
t_c=0.00054;                                                             %cladding thickness
node_f=8; node_c=14;                                                         %node number for fuel side and cladding side. 
node_t=node_f+node_c-1;                                                     %total node number.interface is one between f and c.
node_r=10;          %nodalization R-direction of pin

ddf   =f_o/(node_f-1);                                                      %node distance for fuel
ddc   =t_c/(node_c-1);
ddfh  =H_a/(node_MATRA-1);                                                  %node distance
% r-direction node
n_fr  =0:ddf:f_o;                                                           %z-direction node vector for fuel
n_cr  =f_o+ddc:ddc:f_o+t_c;                                                 %z-direction node vector for clad
n_t   =[n_fr,n_cr];
% r-direction half node
n_frh =[0,ddf/2:ddf:f_o,f_o+ddc/2:ddc:f_o+t_c];                               % r(i+1)/2 node 


DDL=zeros(node_MATRA,node_t);                                            %node volume
DDR=zeros(node_MATRA,node_t);                                            %node volume
KL=zeros(node_MATRA,node_t);                                            %node conductivity
KR=zeros(node_MATRA,node_t);                                            %node conductivity

for i=1:node_MATRA
    for j=1:node_t
        if j<=node_f
            DDL(i,j)=ddf;
        else
            DDL(i,j)=ddc;
        end
    end
    for j=1:node_t
        if j<=node_f-1
            DDR(i,j)=ddf;
        else
            DDR(i,j)=ddc;
        end    
    end
end
V_rf   =zeros(node_MATRA,node_t);                                            %node volume
%% properties
k_uzr   =zeros(node_MATRA,node_t);                                          %FAST pin conductivity
c_uzr   =zeros(node_MATRA,node_t);                                          %heat capacity
den_uzr =zeros(node_MATRA,node_t);                                          %density
%% input value
T_inlet=(390+273);                                                                %Inlet temperature
G=2140;                        
Massflowrate=G*2*A_c_f*217*144;                                            % 217 개의 fuel rod 144 개의 assembly

%% Multi - assembly analysis
Nomalized_O_FA_q=[0.775,0.778,0.755,0.862,0.811,0.864,0.756,0.774,0.775,0.763,0.817,0.762,0.817,0.76,0.907,0.933,0.908,0.911,0.937,0.913,0.87,0.876,0.783,0.761,0.864,0.862,0.923,0.79,0.776,0.91,0.94,0.768,0.776,0.931,0.921,0.758,0.911,0.877,0.865,0.826,0.808,0.873,0.865,0.915,0.915,0.82,0.92,0.918,0.765,0.761,0.773,0.941,0.761,0.76,0.94,0.778,0.778,0.859,0.87,0.778];
Nomalized_I_FA_q=[0.972,0.889,0.981,0.973,1.016,1.024,1.127,1.071,1.129,1.027,1.03,0.989,0.898,0.889,1.07,1.128,1.142,1.217,1.196,1.219,1.154,1.146,1.088,0.987,0.978,1.127,1.197,1.225,1.258,1.276,1.265,1.23,1.209,1.143,1.03,1.022,1.142,1.215,1.254,1.272,1.279,1.276,1.272,1.238,1.158,1.032,1.021,1.123,1.214,1.231,1.272,1.209,1.282,1.228,1.237,1.14,0.992,0.974,1.067,1.199,1.154,1.267,1.082,1.274,1.152,1.208,1.085,0.907,0.894,0.979,1.131,1.032,1.231,0.902,1.233,1.029,1.135,0.986,1.029,1.137,1.139,1.028,0.985,0.987];
Nomalized_q=[Nomalized_O_FA_q,Nomalized_I_FA_q];
Massflux_FA=G*ones(length(Nomalized_q),1);

%% Group
Group_N=20;

Max_G=max(Nomalized_q);
Min_G=min(Nomalized_q);

step_G=(Max_G-Min_G)/(Group_N-1);
Group_d=Min_G-step_G/2:step_G:Max_G+step_G;

Group_q=ones ( 1, length(Nomalized_q));
%  assembly 에 대 해 서 
for i = 1 : length(Nomalized_q)
    for j = 1 : Group_N-1
        if Group_d(j)<=Nomalized_q(i) && Nomalized_q(i) < Group_d(j+1);
            Group_q(i)=mean([Group_d(j),Group_d(j+1)]);
        end
    end
    
    for j = Group_N
        if Group_d(j)<=Nomalized_q(i) && Nomalized_q(i) <= Group_d(j+1);
            Group_q(i)=mean([Group_d(j),Group_d(j+1)]);
        end
    end

end
%Group count
  Group_Q=[];
  
for i = 1 : Group_N
    Group_Q=[Group_Q,Min_G+step_G*(i-1)];                                                               %  Q를 가지는 그룹. N=7일때 step 6 간격으로 가짐.
end
Group_IC=[];
for i = 1 : Group_N-1
    Group_IC=[Group_IC,sum(Min_G+step_G*(i-1)-0.5*step_G<=Group_q(length(Nomalized_O_FA_q)+1:end)&Group_q(length(Nomalized_O_FA_q)+1:end) < Min_G+step_G*(i-1)+0.5*step_G)];          %  Inner Fuel assembly 중에서 Group을 나눔. 1그룹부터 N그룹까지 해당하는 개수를       
end
for i = Group_N
    Group_IC=[Group_IC,sum(Min_G+step_G*(i-1)-0.5*step_G<=Group_q(length(Nomalized_O_FA_q)+1:end)&Group_q(length(Nomalized_O_FA_q)+1:end) <= Min_G+step_G*(i-1)+0.5*step_G)];          %  Inner Fuel assembly 중에서 Group을 나눔. 1그룹부터 N그룹까지 해당하는 개수를       
end
Group_IC=Group_IC*3;                                                                                     % 반응도 삽입하는 갯수를 의미함 3개씩 들어있으므로 곱하기 3 해줘야함.


%% Initial condition

%mass flux[kg/s-m2] %mass flux[kg/s-m2]  2075 는 출구온도를 545 도씨로 만들게 하는 유량이
% T_o should be get from MATRA-LMR (coolant temperature)
node_total=(node_t+1+node_r);                                              %Total node in r-direction

T_uzr = zeros(node_MATRA*node_t,1,length(Group_Q) );
T_coolant=zeros(node_MATRA,1,length(Group_Q));
T_f=zeros(num_h*node_r,1,length(Group_Q));
T_T=zeros(node_total*node_MATRA,1,length(Group_Q));
T_next=zeros(node_total*node_MATRA,1,length(Group_Q));
for AA=1:length(Group_Q)                                               % 모든 어셈블리에 대하여 
T_uzr(:,:,AA) =(390+273)*ones(node_MATRA*node_t,1);                                 %FAST pin temperature region    %steady state condition . same temperature with around the coolant
T_coolant(:,:,AA)=(390+273)*ones(node_MATRA,1);                                    %Initial coolant temperature(MATRA output value)
T_f(:,:,AA)=(390+273)*ones(num_h*node_r,1);                            %FAST pin temperature region
%% TOTAL TEMPERATURE FIELD
T_T(:,:,AA)=zeros(node_total*node_MATRA,1);                                 %Total temperature field(containing FAST and coolant and fuel region)
T_next(:,:,AA)=zeros(node_total*node_MATRA,1);                                 %Total temperature field(containing FAST and coolant and fuel region)

for i=1:node_MATRA
    for j=1:node_total
        if j<=node_t
            T_T(  node_total*(i-1)+j ,1,AA)=T_uzr(j+(i-1)*node_t,1,AA);
        end
        if j==node_t+1
            T_T(  node_total*(i-1)+j ,1,AA)=T_coolant(i,1,AA);
    
        end
        if j>node_t+1
            T_T(  node_total*(i-1)+j ,1,AA)=T_f((j-(node_t+1))+(i-1)*node_r,1,AA);
        end
    end
end
end
%% FAST region temperature transfer time (lumped capacity method)
t_cf=0.0001;        %cladding thickness
num_ah=node_MATRA*(0.4+H_a)/H_a;         %nodalization Z-direction  (Actual core region) up   region
num_dh=node_MATRA*(0.4)/H_a;    %nodalization Z-direction  (Actual core region)  down region    
ddr=D_p/2/(node_r-1);             %node distance
ddh=H_a/(num_h-1);                %node distance
n_r=0:ddr:D_p/2;                  %z-direction node vector
n_rF =[0,ddr/2:ddr:D_p/2];        % r(i+1)/2 node 
n_h=0:ddh:H_a;                    %h-direction node vector
V_r=zeros(num_h,node_r);                                       %node volume


%steady state condition . same temperature with around the coolant
%% size setting for calculating speed



plota=zeros(5000,1);
plotv=zeros(5000,1);
plots=zeros(5000,1);
plotb=zeros(5000,1);
plotg=zeros(5000,1);
plotd=zeros(5000,1);
plotq=zeros(5000,1);              
plotp=zeros(5000,1);
plotpd=zeros(5000,1);
plot_position_start=ones(5000,1)*zi_transient_F(1);
plot_position_end=ones(5000,1)*zi_transient_F(end);
plot_absorber_h=ones(5000,1)*ab_transient_F;
plot_reactivity=ones(5000,1);
plot_feedback=zeros(5000,1);
plot_inserted_r=zeros(5000,1);
plot_Tem_f=zeros(5000,1);
plot_Tem_c=zeros(5000,1);
plot_Tem_s=zeros(5000,1);
plot_Tem_cm=zeros(5000,1);
plot_Tem_fm=zeros(5000,1);
plott=zeros(5000,1);
plot_den_b=zeros(5000,Group_N);                                                  % density at the bottom
plot_den_t=zeros(5000,Group_N);                                                  % density at the top
plot_Tcoolant=zeros(5000,3,Group_N);                                               % Coolant  temperature Z = Ha (active core height), Z=Ha/2, Z= 0
plot_Tpincoolant=zeros(5000,3,Group_N);                                            % pin Cooalnt temperature the pin Z = Ha (active core height), Z=Ha/2, Z= 0
plot_Tpin=zeros(5000,3,Group_N);                                                   % pin  temperature Z = Ha (active core height), Z=Ha/2, Z= 0
plot_Tfuel=zeros(5000,3,Group_N);                                                  % Fuel_wall  temperature Z = Ha (active core height), Z=Ha/2, Z= 0
plot_Tcladding=zeros(5000,3,Group_N);                                              % cladding_wall temperature the pin Z = Ha (active core height), Z=Ha/2, Z= 0
plot_xe = zeros(5000,Group_N);



c1=zeros(1,num_gap-1);
Shear_Drag=zeros(num_h-1,1);                                          %size setting
Fb=zeros(num_h-1,1);


k_coolant=zeros(node_MATRA ,1);
den_coolant=zeros(node_MATRA ,1,Group_N);
cp_coolant=  zeros(node_MATRA ,1);
mu_coolant=zeros(node_MATRA ,1);
Re=zeros(node_MATRA ,1,Group_N);
Pr=zeros(node_MATRA ,1);
H_c=zeros(node_MATRA ,1);


Plot_T=zeros(node_total,node_MATRA,Group_N);
Plot_TA=zeros(1+node_t,node_MATRA,Group_N);
Plot_TB=zeros(node_total-node_t-1,node_MATRA,Group_N);
Plot_TC=zeros(node_total,node_MATRA,Group_N);
k_f=zeros(num_h,node_r);                                       %FAST pin conductivity
c_f=zeros(num_h,node_r);                                       %heat capacity
den_f=zeros(num_h,node_r);                                     %density



A_l=zeros(num_h,node_r);    A_r=zeros(num_h,node_r);          %node right side area      %node left side area
A_u=zeros(num_h,node_r);    A_d=zeros(num_h,node_r);          %node up side area      %node down side area
Ldt=0;

A_lf=zeros(node_MATRA,node_t);    A_rf=zeros(node_MATRA,node_t);             %node right side area      %node left side area
A_uf=zeros(node_MATRA,node_t);    A_df=zeros(node_MATRA,node_t);             %node up side area      %node down side area
q_m=zeros(node_MATRA,node_t);                                              %matrics of heat
M = zeros((node_r+node_t+1)*node_MATRA,(node_r+node_t+1)*node_MATRA); 
MH= zeros((node_r+node_t+1)*node_MATRA,1);                 % heat generation vector
step=0;
syms x;





error=1;
while error>1e-8

    step=step+1;
    time=dt*step;
    
    
for AA = 1: length(Group_Q)

for j=1:node_t                                                              %for overall region including fuel and clad
    for i=1:node_MATRA
        % volume & area
        if j < node_f                                                      % Fuel 영역에서
            V_rf(i,j)=2*pi*ddfh*ddf*n_t(j);                                % unit volume
            A_rf(i,j)=2*pi*(n_t(j)+ddf/2)*ddfh;                            % node right side area
            A_lf(i,j)=2*pi*(n_t(j)-ddf/2)*ddfh;                            % node left side area
            A_uf(i,j)=2*pi*ddf*n_t(j);                                     % node up side area
            A_df(i,j)=2*pi*ddf*n_t(j);                                     % node up side area
            if j==1                                                        % 왼쪽 단열 (중심 부에서)
                V_rf(i,j)=pi/4*ddf^2*ddfh;
                A_lf(i,j)=0;                         
                A_uf(i,j)=pi/4*ddf^2;
                A_df(i,j)=A_uf(i,j);
            end

            if i==1                                                             %맨 윗면에 대해서 면적 및 부피가 기존것의 반. dh가 절반이기 때문
                V_rf(i,j)=V_rf(i,j)/2;
                A_rf(i,j)=A_rf(i,j)/2;
                A_lf(i,j)=A_lf(i,j)/2;
                A_uf(i,j)=0;                                                     %위쪽방향으로 열전달 x 
            end
            if i==node_MATRA                                                    %맨 아랫면에 대해서 면적 및 부피가 기존의 반. dh가 절반이기 때문
                V_rf(i,j)=V_rf(i,j)/2;
                A_rf(i,j)=A_rf(i,j)/2;
                A_lf(i,j)=A_lf(i,j)/2;
                A_df(i,j)=0;                                                     %아래방향으로 열전달 x
            end 
        end


        if j > node_f                                                      %Clad 영역에서
            V_rf(i,j)=2*pi*ddfh*ddc*n_t(j);                                %unit volume
            A_rf(i,j)=2*pi*(n_t(j)+ddc/2)*ddfh;                            %node right side area
            A_lf(i,j)=2*pi*(n_t(j)-ddc/2)*ddfh;                            %node left side area
            A_uf(i,j)=2*pi*ddc*n_t(j);                                     %node up side area
            A_df(i,j)=2*pi*ddc*n_t(j);                                     %node up side area

            if j==node_t                                                   %맨 오른쪽 경계면에서 ( 반쪼가리)
                V_rf(i,j)=pi*ddfh*ddc*(n_t(j)-ddc/4);
                A_rf(i,j)=2*pi*(n_t(j))*ddfh;
                A_uf(i,j)=pi*ddc*(n_t(j)-ddc/4);
                A_df(i,j)=A_uf(i,j);
            end
            if i==1                                                        %맨 윗면에 대해서 면적 및 부피가 기존것의 반. dh가 절반이기 때문
                V_rf(i,j)=V_rf(i,j)/2;
                A_rf(i,j)=A_rf(i,j)/2;
                A_lf(i,j)=A_lf(i,j)/2;
                A_uf(i,j)=0;                                               %위쪽방향으로 열전달 x 
            end
            if i==node_MATRA                                               %맨 아랫면에 대해서 면적 및 부피가 기존의 반. dh가 절반이기 때문
                V_rf(i,j)=V_rf(i,j)/2;
                A_rf(i,j)=A_rf(i,j)/2;
                A_lf(i,j)=A_lf(i,j)/2;
                A_df(i,j)=0;                                               %아래방향으로 열전달 x
            end 
        end


        if j == node_f                                                         %경계면에서 Fuel and Cladding interface
            V_ff = pi*ddfh*ddf*(n_t(j)-ddf/4) ;                               %volume fraction of fuel
            V_fc = pi*ddfh*ddc*(n_t(j)+ddc/4) ;                                %volume fraction of cla
            
            V_rf(i,j)=V_ff+V_fc;    %unit volume
            A_rf(i,j)=2*pi*(n_t(j)+ddc/2)*ddfh;                                 %node right side area
            A_lf(i,j)=2*pi*(n_t(j)-ddf/2)*ddfh;                                 %node left side area
            A_uf(i,j)=pi*ddf*(n_t(j)-ddf/4)+pi*ddc*(n_t(j)+ddc/4);              %node up side area
            A_df(i,j)=pi*ddf*(n_t(j)-ddf/4)+pi*ddc*(n_t(j)+ddc/4);              %node up side area
            if i==1                                                             %맨 윗면에 대해서 면적 및 부피가 기존것의 반. dh가 절반이기 때문
                V_rf(i,j)=V_rf(i,j)/2;
                A_rf(i,j)=A_rf(i,j)/2;
                A_lf(i,j)=A_lf(i,j)/2;
                A_uf(i,j)=0;                                                     %위쪽방향으로 열전달 x 
            end
            if i==node_MATRA                                                    %맨 아랫면에 대해서 면적 및 부피가 기존의 반. dh가 절반이기 때문
                V_rf(i,j)=V_rf(i,j)/2;
                A_rf(i,j)=A_rf(i,j)/2;
                A_lf(i,j)=A_lf(i,j)/2;
                A_df(i,j)=0;                                                     %아래방향으로 열전달 x
            end 
        end        

        %%properties
        %conductivity, density, capacity
        if j<=node_f                                      %Uzr region
            k_uzr(i,j)=S(T_uzr((i-1)*node_t+j,1,AA));     %[W/m.K]
            den_uzr(i,j)=12375.864; %density matrix
            c_uzr(i,j)=121.336;                                        %capacity matrix
        else                                            %HT-9 region
            k_uzr(i,j)=17.622+2.428*10^(-2)*(T_uzr((i-1)*node_t+j,1,AA))-1.696*10^(-5)*(T_uzr((i-1)*node_t+j,1,AA))^2;%[W/m.K]
            if  T_uzr((i-1)*node_t+j,1,AA)<800
                c_uzr(i,j)=561;%(T_uzr((i-1)*node_t+j,1)-500)/6+500;
            else
                c_uzr(i,j)=561;%(T_uzr((i-1)*node_t+j,1)-800)/6+550;
            end
             den_uzr(i,j)=7750;% 7700~7800 ferritic S/S 
        end
       
        %properties at the interface follows volume fraction.
        if j==node_f
            k_uzr(i,j)=(S(T_uzr((i-1)*node_t+j,1,AA)) + 17.622+2.428*10^(-2)*(T_uzr((i-1)*node_t+j,1,AA))-1.696*10^(-5)*(T_uzr((i-1)*node_t+j,1,AA))^2 )/2;      %[W/m.K]

        end
        
        
    end
end

%% 좌변/우변에 곱해지는 특성치들을 고려하기 위해서   KL = 왼쪽 노드, KR= 오른쪽 노드
for i=1:node_MATRA
    for j=1:node_t
        if j<=node_f
            KL(i,j)=S(T_uzr((i-1)*node_t+j,1,AA));                                %[W/m.K] Fuel side
        else
            KL(i,j)=17.622+2.428*10^(-2)*(T_uzr((i-1)*node_t+j,1,AA))-1.696*10^(-5)*(T_uzr((i-1)*node_t+j,1,AA))^2;   %[W/m.K] cladding side
        end
    end
    for j=1:node_t
        if j<=node_f-1
            KR(i,j)=S(T_uzr((i-1)*node_t+j,1,AA));                                %[W/m.K] Fuel side
        else
            KR(i,j)=17.622+2.428*10^(-2)*(T_uzr((i-1)*node_t+j,1,AA))-1.696*10^(-5)*(T_uzr((i-1)*node_t+j,1,AA))^2;   %[W/m.K] cladding side
        end    
    end
end

%% 조화 평균 conductivity
% KL_m=KL; % 변수 선언
% KR_m=KR; % 변수 선언
% % Centerline에 대해서
% KL_m(1:end,1)=0;                                                           % 중앙에서 왼쪽으로 열전달 x
% KR_m(1:end,1)=k_uzr(1:end,1);                                              % 오른쪽으로 열 전달시 conductivity는 centerline것을 이용. ( 추후 변경 가능 k i+1/2 로)
% % 조화평균 사용                                                             %
% % km(2*k1k2/(k1+k2)
% for i=2:length(KL(1,1:end))
%     k1=KL(1:end,i-1);
%     k2=KL(1:end,i);
%     KL_m(1:end,i)=2*k1.*k2./(k1+k2);                                       % 좌변 노드에 대한 물성치
%     
%     k1=KR(1:end,i-1);
%     k2=KR(1:end,i);
%     KR_m(1:end,i)=2*k1.*k2./(k1+k2);                                       % 우변 노드에 대한 물성치
% end
%Conductivity
k_uzr_m=zeros(length(k_uzr(1:end,1)),length(k_uzr(1,1:end))-1);            %매트릭스 생성
for i=2:length(KL(1,1:end))
    k1=k_uzr(1:end,i-1);
    k2=k_uzr(1:end,i);
    k_uzr_m(1:end,i-1)=2*k1.*k2./(k1+k2); 
end
%중심에서 열전도계수
k_uzr_m=[k_uzr(1:end,1),k_uzr_m];

%heat flux
A_q=A_df(2,1:node_f);  A_q(1,node_f)=pi*ddf*(n_t(node_f)-ddf/4);                           %heat flux에 곱해질 수
n_h=[H_a/2,(H_a/2-ddfh/2):-ddfh:(-H_a/2+ddfh/2),-H_a/2];                    %heat flux에 적분될 변수들
% q_m = heat generation
for i=1:node_MATRA
    for j=1:node_f                                                         %for overall region including fuel and clad
        q_m(i,j)=A_q(j)/pi*con_q*H_a*(sin(pi/H_a*n_h(i))-sin(pi/H_a*n_h(i+1))) * (Group_Q(AA)) ;  % heat generation node
    end
end
% q_mv=volumetric heat generation
q_mv=q_m./V_rf;


%%Active core region coolant side
for i=1:node_MATRA
%% liquid sodium properties
    k_coolant(i,1)=92.951-5.8087e-2*(T_coolant(i,1,AA)-273.15)+11.7274e-6*(T_coolant(i,1,AA)-273.15)^2;                                       %k matrix
    den_coolant(i,1,AA)=950.0483-0.2298*(T_coolant(i,1,AA)-273.15)-14.6045e-6*(T_coolant(i,1,AA)-273.15)^2+5.6377e-9*(T_coolant(i,1,AA)-273.15)^3; %density matrix
    cp_coolant(i,1)=  1436.715-0.5806*(T_coolant(i,1,AA)-273.15)+4.6273e-4*(T_coolant(i,1,AA)-273.15)^2; %          Coolant_cp(T_coolant(i,1)-273.15); %                           %capacity matrix
    mu_coolant(i,1)=10^(-2.4892+220.65/T_coolant(i,1,AA)-0.4295*log10(T_coolant(i,1,AA)));                                                    %mu matrix   
    Re(i,1,AA)=G*D_fh/mu_coolant(i,1);                                                                                              %Reynolds number
    Pr(i,1)=cp_coolant(i,1)*mu_coolant(i,1)/k_coolant(i,1);                                                                                         %Plantle number
    H_c(i,1)=k_coolant(i,1)/D_fh*0.023*Re(i,1)^0.8*Pr(i,1)^0.4;                                                                   %Dittus-Boelter equation. heat transfer coefficient

end
for j=1:node_r
    for i=1:num_h
        % volume & area
        V_r(i,j)=2*pi*ddh*ddr*n_r(j);           %unit volume
        A_r(i,j)=2*pi*(n_r(j)+ddr/2)*ddh;       %node right side area
        A_l(i,j)=2*pi*(n_r(j)-ddr/2)*ddh;       %node left side area
        A_u(i,j)=2*pi*ddr*n_r(j);               %node up side area
        A_d(i,j)=2*pi*ddr*n_r(j);               %node up side area
        if j==1
            V_r(i,j)=pi/4*ddr^2*ddh;
            A_l(i,j)=0;                         %왼쪽 단열
            A_u(i,j)=pi/4*ddr^2;
            A_d(i,j)=A_u(i,j);
        end
        if j==node_r
            V_r(i,j)=pi*ddh*ddr*(n_r(j)-ddr/4);
            A_r(i,j)=2*pi*(n_r(j))*ddh;
            A_u(i,j)=pi*ddr*(n_r(j)-ddr/4);
            A_d(i,j)=A_u(i,j);
        end
        if i==1
            V_r(i,j)=V_r(i,j)/2;
            A_r(i,j)=A_r(i,j)/2;
            A_l(i,j)=A_l(i,j)/2;
            A_u(i,j)=0;
        end
        if i==num_h
            V_r(i,j)=V_r(i,j)/2;
            A_r(i,j)=A_r(i,j)/2;
            A_l(i,j)=A_l(i,j)/2;
            A_d(i,j)=0;
        end 
         % properties
         %conductivity
        if j<=floor((D_p/2-t_cf)/ddr)                    %fluid region
            k_f(i,j)=92.951-5.8087e-2*(T_f((i-1)*node_r+j,1,AA)-273.15)+11.7274e-6*(T_f((i-1)*node_r+j,1,AA)-273.15)^2;     %[W/m.K]
            den_f(i,j)=950.0483-0.2298*(T_f((i-1)*node_r+j,1,AA)-273.15)-14.6045e-6*(T_f((i-1)*node_r+j,1,AA)-273.15)^2+5.6377e-9*(T_f((i-1)*node_r+j,1,AA)-273.15)^3; %density matrix
            c_f(i,j)=1436.715-0.5806*(T_f((i-1)*node_r+j,1,AA)-273.15)+4.6273e-4*(T_f((i-1)*node_r+j,1,AA)-273.15)^2;    %Coolant_cp((T_f((i-1)*node_r+j,1)-273.15));% 
            
            %                                   %capacity matrix
        else                                            %HT-9 region
            k_f(i,j)=17.622+2.428*10^(-2)*(T_f((i-1)*node_r+j,1,AA))-1.696*10^(-5)*(T_f((i-1)*node_r+j,1,AA))^2;%[W/m.K]
            if T_f((i-1)*node_r+j,1)<800
                c_f(i,j)=(T_f((i-1)*node_r+j,1,AA)-500)/6+500;
            else
                c_f(i,j)=(T_f((i-1)*node_r+j,1,AA)-800)/6+550;
            end
             den_f(i,j)=7750;% 7700~7800 ferritic S/S 
        end
    end
end

%% Fuel and cladding Matrix
mcf =den_uzr.*V_rf.*c_uzr./dt;                                             %Matrix coefficient (alpha)
% boundary condition at the interface
density_fuelside= 12375.864*ones(node_MATRA,1);     %[kg/m3]
density_cladside= 7750*ones(node_MATRA,1);% 
volume_fuelside = pi*ddfh*ddf*(n_t(node_f)-ddf/4)*ones(node_MATRA,1) ;                               %volume fraction of fuel
volume_cladside = pi*ddfh*ddc*(n_t(node_f)+ddc/4)*ones(node_MATRA,1) ;                               %volume fraction of cla
volume_fuelside(1,1)=volume_fuelside(1,1)/2; volume_fuelside(end,1)=volume_fuelside(end,1)/2;        %end point is half
volume_cladside(1,1)=volume_cladside(1,1)/2; volume_cladside(end,1)=volume_cladside(end,1)/2;
capacity_fuelside=121.336;
capacity_cladside=561;

%pcp/dt*V
mcf(1:end,node_f)=(density_fuelside.*volume_fuelside*capacity_fuelside./dt)+(density_cladside.*volume_cladside*capacity_cladside./dt);
%pcp/dt
vcf=mcf./V_rf;

%% FAST region matrix
mc =den_f.*V_r.*c_f./dt;            
vc =mc./V_r;


%% matrix 생성
% Interior node 에 대해서 먼저 실행
NN=node_t;
MM=node_MATRA;                                % z 방향 노드
FF=node_r;                                    % FAST r 방향 노드
CN=1;                                         % 쿨런트 r 방향 노드
TT=NN+FF+CN;                                  % 총 노드




%%  Interior node
for i=1*TT+1:(MM-1)*TT-1
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end

    % FUEL ROD REGION
    if (node_j>1) && (node_j<NN)
        %heat generation
        MH(i,1)=q_mv(node_i,node_j)/vcf(node_i,node_j);
        % Tob ( conductivity 는 중앙 노드꺼씀)
        TC=-k_uzr(node_i,node_j)/ddfh^2/vcf(node_i,node_j);                        
        % BOT ( conductivity 는 중앙 노드꺼씀)
        BC=-k_uzr(node_i,node_j)/ddfh^2/vcf(node_i,node_j);                         
        % Left (조화평균값)
        LC=-n_frh(node_j)*k_uzr_m(node_i,node_j)    /n_t(node_j)/(n_frh(node_j+1)-n_frh(node_j))/vcf(node_i,node_j)/(n_t(node_j)  -n_t(node_j-1));
        % Right (conductivity 조화 평균)
        RC=-n_frh(node_j+1)*k_uzr_m(node_i,node_j+1)/n_t(node_j)/(n_frh(node_j+1)-n_frh(node_j))/vcf(node_i,node_j)/(n_t(node_j+1)-n_t(node_j));
        % Center
        CC=-(TC+BC+LC+RC)+1;

        M(i,i-TT)=TC;
        M(i,i+1)=RC;
        M(i,i)=CC;
        M(i,i-1)=LC;
        M(i,i+TT)=BC;
    
    end
    
        % FAST region
    if (node_j>NN+CN+1) && (node_j<TT)
        node_j=node_j-(NN+CN);                                             %인덱스 조정
        % Tob ( conductivity 는 중앙 노드꺼씀)
        TC=-k_f(node_i,node_j)/ddh^2/vc(node_i,node_j);                        
        % BOT ( conductivity 는 중앙 노드꺼씀)
        BC=-k_f(node_i,node_j)/ddh^2/vc(node_i,node_j);                         
        % Left (조화평균값)
        LC=-n_rF(node_j)*k_f(node_i,node_j)    /n_r(node_j)/(n_rF(node_j+1)-n_rF(node_j))/vc(node_i,node_j)/(n_r(node_j)  -n_r(node_j-1));
        % Right (conductivity 조화 평균)
        RC=-n_rF(node_j+1)*k_f(node_i,node_j)/n_r(node_j)/(n_rF(node_j+1)-n_rF(node_j))/vc(node_i,node_j)/(n_r(node_j+1)-n_r(node_j));
        % Center
        CC=-(TC+BC+LC+RC)+1;

        M(i,i-TT)=TC;
        M(i,i+1)=RC;
        M(i,i)=CC;
        M(i,i-1)=LC;
        M(i,i+TT)=BC;
    
    end
    
    
end
%%  Fuel centerline node
for i=1:TT:TT*MM
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    %heat generation
    MH(i,1)=q_mv(node_i,node_j)/vcf(node_i,node_j);

    % Right (conductivity 조화 평균)
    RC=-4*k_uzr_m(node_i,node_j)/vcf(node_i,node_j)/(n_t(node_j+1)-n_t(node_j))^2;
    % Center
    CC=-(RC)+1;
    
    M(i,i+1)=RC;
    M(i,i)=CC;
end
%% FAST centerline node
for i=NN+CN+1:TT:TT*MM
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    node_j=node_j-(NN+CN);                                                 %인덱스 조정
    %  (conductivity 일반 값 사용)
    RC=-4*k_f(node_i,node_j)/vc(node_i,node_j)/(n_r(node_j+1)-n_r(node_j))^2;
    % Center
    CC=-(RC)+1;
    
    M(i,i+1)=RC;
    M(i,i)=CC;
    
end

%%  Fuel upper node
for i=2:NN-1
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    %heat generation
    MH(i,1)=q_mv(node_i,node_j)/vcf(node_i,node_j);

    % BOT ( conductivity 는 중앙 노드꺼씀)
    BC=-k_uzr(node_i,node_j)/ddfh^2/vcf(node_i,node_j);                         
    % Left (조화평균값)
    LC=-n_frh(node_j)*k_uzr_m(node_i,node_j)/n_t(node_j)/(n_frh(node_j+1)-n_frh(node_j))/vcf(node_i,node_j)/(n_t(node_j)-n_t(node_j-1));
    % Right (conductivity 조화 평균)
    RC=-n_frh(node_j+1)*k_uzr_m(node_i,node_j+1)/n_t(node_j)/(n_frh(node_j+1)-n_frh(node_j))/vcf(node_i,node_j)/(n_t(node_j+1)-n_t(node_j));
    % Center
    CC=-(BC+LC+RC)+1;
  

    M(i,i+1)=RC;
    M(i,i)=CC;
    M(i,i-1)=LC;
    M(i,i+TT)=BC;
    
end
%%  FAST upper node
for i=NN+CN+2:TT-1
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    node_j=node_j-(NN+CN);                                             %인덱스 조정
    % BOT ( conductivity 는 중앙 노드꺼씀)
    BC=-k_f(node_i,node_j)/ddh^2/vc(node_i,node_j);                         
    % Left (조화평균값)
    LC=-n_rF(node_j)*k_f(node_i,node_j)/n_r(node_j)/(n_rF(node_j+1)-n_rF(node_j))/vc(node_i,node_j)/(n_r(node_j)-n_r(node_j-1));
    % Right (conductivity 조화 평균)
    RC=-n_rF(node_j+1)*k_f(node_i,node_j)/n_r(node_j)/(n_rF(node_j+1)-n_rF(node_j))/vc(node_i,node_j)/(n_r(node_j+1)-n_r(node_j));
    % Center
    CC=-(BC+LC+RC)+1;
  

    M(i,i+1)=RC;
    M(i,i)=CC;
    M(i,i-1)=LC;
    M(i,i+TT)=BC;
    
end

%%  Fuel bottom node
for i=(MM-1)*(TT)+2:MM*TT-(FF+2)
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    %heat generation
    MH(i,1)=q_mv(node_i,node_j)/vcf(node_i,node_j);
    % Tob ( conductivity 는 중앙 노드꺼씀)
    TC=-k_uzr(node_i,node_j)/ddfh^2/vcf(node_i,node_j);   
    % Left (조화평균값)
    LC=-n_frh(node_j)*k_uzr_m(node_i,node_j)    /n_t(node_j)/(n_frh(node_j+1)-n_frh(node_j))/vcf(node_i,node_j)/(n_t(node_j)  -n_t(node_j-1));
    % Right (conductivity 조화 평균)
    RC=-n_frh(node_j+1)*k_uzr_m(node_i,node_j+1)/n_t(node_j)/(n_frh(node_j+1)-n_frh(node_j))/vcf(node_i,node_j)/(n_t(node_j+1)-n_t(node_j));
    % Center
    CC=-(TC+LC+RC)+1;
  
    M(i,i-TT)=TC;
    M(i,i+1)=RC;
    M(i,i)=CC;
    M(i,i-1)=LC;
end
%%  FAST bottom node
for i=(MM-1)*(TT)+(NN+CN)+2:MM*TT-(1)
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    node_j=node_j-(NN+CN);                                             %인덱스 조정

    % Tob ( conductivity 는 중앙 노드꺼씀)
    TC=-k_f(node_i,node_j)/ddh^2/vc(node_i,node_j);   
    % Left (조화평균값)
    LC=-n_rF(node_j)*k_f(node_i,node_j)    /n_r(node_j)/(n_rF(node_j+1)-n_rF(node_j))/vc(node_i,node_j)/(n_r(node_j)  -n_r(node_j-1));
    % Right (conductivity 조화 평균)
    RC=-n_rF(node_j+1)*k_f(node_i,node_j)/n_r(node_j)/(n_rF(node_j+1)-n_rF(node_j))/vc(node_i,node_j)/(n_r(node_j+1)-n_r(node_j));
    % Center
    CC=-(TC+LC+RC)+1;
  
    M(i,i-TT)=TC;
    M(i,i+1)=RC;
    M(i,i)=CC;
    M(i,i-1)=LC;
end

%%  Fuel rod outter node
for i=NN:TT:TT*MM
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    %heat generation
    MH(i,1)=q_mv(node_i,node_j)/vcf(node_i,node_j);

    % Tob ( conductivity 는 중앙 노드꺼씀)
    TC=-k_uzr(node_i,node_j)/ddfh^2/vcf(node_i,node_j);                        
    % BOT ( conductivity 는 중앙 노드꺼씀)
    BC=-k_uzr(node_i,node_j)/ddfh^2/vcf(node_i,node_j);                         
    % Left (조화평균값)
    LC=-k_uzr_m(node_i,node_j)*A_lf(node_i,node_j)/(n_t(node_j)-n_t(node_j-1))/mcf(node_i,node_j);
    % Right (conductivity 조화 평균)
    RC=-H_c(node_i,1)*A_rf(node_i,node_j)/mcf(node_i,node_j);
    CC=-(TC+BC+LC+RC)+1;

        if node_i==1
            CC=-(BC+LC+RC)+1;                                                                   % 맨 위쪽 열전달 x
            M(i,i+1)=RC;
            M(i,i)=CC;
            M(i,i-1)=LC;
            M(i,i+TT)=BC;
        end
        if (1<node_i) && (node_i<MM)
            CC=-(TC+BC+LC+RC)+1;
            M(i,i-TT)=TC;
            M(i,i+1)=RC;
            M(i,i)=CC;
            M(i,i-1)=LC;
            M(i,i+TT)=BC;
        end
        if node_i==MM                                                        % 맨 아래쪽 열전달 x
            CC=-(TC+LC+RC)+1;
            M(i,i-TT)=TC;
            M(i,i+1)=RC;
            M(i,i)=CC;
            M(i,i-1)=LC;
        end 

    

end

%%  FAST pin outter node
for i=TT:TT:TT*MM
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end

    node_j=node_j-(NN+CN);                                             %인덱스 조정

    % Tob ( conductivity 는 중앙 노드꺼씀)
    TC=-k_f(node_i,node_j)/ddh^2/vc(node_i,node_j);                        
    % BOT ( conductivity 는 중앙 노드꺼씀)
    BC=-k_f(node_i,node_j)/ddh^2/vc(node_i,node_j);                         
    % Left (조화평균값)
    LC=-k_f(node_i,node_j)*A_l(node_i,node_j)/(n_r(node_j)-n_r(node_j-1))/mc(node_i,node_j);
    % Right (conductivity 조화 평균)
    RC=-H_c(node_i,1)*A_r(node_i,node_j)/mc(node_i,node_j);
    CC=-(TC+BC+LC+RC)+1;

        if node_i==1
            CC=-(BC+LC+RC)+1;                                                                   % 맨 위쪽 열전달 x
            M(i,i-FF)=RC;
            M(i,i)=CC;
            M(i,i-1)=LC;
            M(i,i+TT)=BC;
        end
        if (1<node_i) && (node_i<MM)
            CC=-(TC+BC+LC+RC)+1;
            M(i,i-TT)=TC;
            M(i,i-FF)=RC;
            M(i,i)=CC;
            M(i,i-1)=LC;
            M(i,i+TT)=BC;
        end
        if node_i==MM                                                        % 맨 아래쪽 열전달 x
            CC=-(TC+LC+RC)+1;
            M(i,i-TT)=TC;
            M(i,i-FF)=RC;
            M(i,i)=CC;
            M(i,i-1)=LC;
        end 

end

%% Coolant hyperbolic equation
for nodei=NN+CN:TT:TT*MM
    j=rem(nodei,TT);                  % node j index
    i=fix(nodei/TT)+1;                % node i index
    if j==0
        j=TT;
        i=i-1;
    end
%% Hyperboilic equation (Coolant matrix)
a1=Massflux_FA(AA)*dt/den_coolant(i,1,AA)/(n_h(1,i)-n_h(1,i+1));                            %Coolant in and out energy transfer
b1=1/2*A_rf(i,end)*dt/den_coolant(i,1,AA)/(A_c_f)/(n_h(1,i)-n_h(1,i+1))/cp_coolant(i,1)*H_c(i,1);  %heat source-coolant interaction ( 1/3 = 120/360)
c1=1/6*A_rf(i,end)*dt/den_coolant(i,1,AA)/(A_c_f)/(n_h(1,i)-n_h(1,i+1))/cp_coolant(i,1)*H_c(i,1);  %fast-coolant interaction (1/6 = 60/360 )
M(nodei,nodei)=1+a1+b1+c1;                                                 % - center matrix variable
if i<MM
M(nodei,nodei+TT)=-a1;                                                     % - previous coolant matrix variable
else
MH(nodei,1)=a1*T_inlet;                                                   %inlet boundary
end
M(nodei,nodei-1)=-b1;                                                      % - cladding - coolant
M(nodei,nodei+FF)=-c1;                                                     % - fast - coolant

end

T_next(:,:,AA)=M\T_T(:,:,AA)+M\MH;

end




error=abs(norm(T_T(:,:,1))-norm(T_next(:,:,1)))/norm(T_T(:,:,1));
for AA=1:length(Group_Q)

%% Error term

for i=1:node_MATRA
    for j=1:node_total
        if j<=node_t
            T_uzr(j+(i-1)*node_t,1,AA)=T_T(  node_total*(i-1)+j ,1,AA);
        end
        if j==node_t+1
            T_coolant(i,1,AA)=T_T(  node_total*(i-1)+j ,1,AA);
    
        end
        if j>node_t+1
            T_f((j-(node_t+1))+(i-1)*node_r,1,AA)=T_T(  node_total*(i-1)+j ,1,AA);
        end
    end
end
T_T(:,:,AA)=T_next(:,:,AA);


end
'move to next time step';
%    figure (1)
for AA = 1 : length(Group_Q)
   Plot_T(:,:,AA)=reshape(T_next(:,:,AA),node_total,node_MATRA);
   Plot_TA(:,:,AA)=Plot_T(1:node_t+1,:,AA);
   Plot_TB(:,:,AA)=flipud(Plot_T(node_t+2:end,:,AA));
   Plot_TC(:,:,AA)=[Plot_TA(:,:,AA);Plot_TB(:,:,AA)];

%  AAA=sum((Plot_TC(node_t,1:end)-Plot_TC(node_t+1,1:end)).*(reshape(H_c,1,node_MATRA)).*reshape(A_r(1:end,end),1,node_MATRA))
%  BBB=G*A_c_f*2*cp_coolant(1,1)*(Plot_TC(node_t+1,1)-Plot_TC(node_t+1,end))
end
%surf(Plot_TC(:,:,AA))
%plot(Plot_TC(:,1,AA))
%drawnow;
end
Plot_TC
node_total
node_MATRA
'Assembly temperature done'
%% Pressure loss in assembly
temporary_integral=0;
Buoyancy_c = zeros(1,Group_N);
for AA=1:length(Group_Q)
Buoyancy_c (:,AA)=sum(den_coolant(1:end-1,:,AA)*(-g)*ddh);               % bouyancy term
% contraction and expansion coefficient
K_expansion = (1-(pi*D_p^2/4/2)/(3^0.5/4*(D_p*PD)^2))^2;
K_contraction = 0.42*(1-(pi*D_p^2/4/2)/(3^0.5/4*(D_p*PD)^2));
%CTD model 

F = (0.8063-0.9022*log10(H_a/(D_p+D_w))+0.3526*(log10(H_a/(D_p+D_w)))^2)*(PD)^9.7*(H_a/(D_p+D_w))^(1.78-2*(PD));                        %coefficient for friction factor
f = mean( F./Re(:,:,AA).^0.18);


K(:,AA)= 0.5*K_expansion/den_coolant(1,:,AA)+0.5*K_contraction/den_coolant(end,:,AA)+f*(H_a)/D_fh*0.5/mean(den_coolant(1,:,AA));

temporary =(( x-Buoyancy_c(:,AA))/(K(:,AA)))^0.5;
temporary_integral=temporary_integral+temporary;

eachpressuredrop (AA,1)=Buoyancy_c (:,AA)+0.5*K_expansion/den_coolant(1,:,AA)*G^2++0.5*K_contraction/den_coolant(end,:,AA)*G^2+f*(H_a)/D_fh*0.5/mean(den_coolant(1,:,AA))*G^2;

end
% problem solve  ( assembly에 흐르는 유량 ) 
% equation =temporary_integral==Massflowrate/A_c_f/2/217/144;

%pressure drop
% PressureDrop=double(solve(equation));

% Massflux_FA=((PressureDrop-Buoyancy_c (:,:))./K(:,:)).^0.5;
% 'mass flux distribution done'
%% Steady-state complete.

%% for speed
F_buoyancy=zeros(1,Group_N); F_drag=zeros(1,Group_N); 
F_pressure=zeros(1,Group_N);
Q_error=0;
Pe_dz=zeros(1,Group_N);
a=zeros(1,Group_N);
v=zeros(1,Group_N);
s=zeros(1,Group_N);
trigger=zeros(1,Group_N);
q_t=zeros(1,Group_N);
v_i=zeros(num_gap,1,Group_N);                                               
v_i_p=zeros(num_gap,1,Group_N);         
con_q=zeros(1,Group_N);
dP=zeros(1,Group_N);



%% Initial Condition setting
for AA = 1 : length(Group_Q)
zi_transient_F(:,:,AA)=zi_F+(H_p-H);                                        %actual position of FAST per height    (vector )
ab_transient_F=H_p-h_ab;                                                   %Absorber height from the bottom    ( scalar ) 
Fb=zeros(node_MATRA-1,AA);

    

Interpx=H_a:-ddfh:0;                                %인터폴레이션 x 축
Interpy=Plot_TC(end-round(R/ddr),1:end,AA);            %FAST 벽면에 해당하는 온도    
for i=1:node_MATRA                              %interpolation을 통해 온도 매칭.
    T_c_steady(i,1,AA)=interp1(Interpx,Interpy,(i-1)*dz);
    if (i-1)*dz>H_a
        T_c_steady(i,1,AA)=T_c_steady(i-1,1,AA);
    end       
end
for i=1:node_MATRA
%liquid sodium properties
    k_c(i,1)=92.951-5.8087e-2*(T_c_steady(i,1,AA)-273.15)+11.7274e-6*(T_c_steady(i,1,AA)-273.15)^2;                                       %k matrix
    den_c(i,1)=950.0483-0.2298*(T_c_steady(i,1,AA)-273.15)-14.6045e-6*(T_c_steady(i,1,AA)-273.15)^2+5.6377e-9*(T_c_steady(i,1,AA)-273.15)^3; %density matrix
    cp_c(i,1)=1436.715-0.5806*(T_c_steady(i,1)-273.15)+4.6273e-4*(T_c_steady(i,1)-273.15)^2;                                        %capacity matrix
    mu_c(i,1)=10^(-2.4892+220.65/T_c_steady(i,1,AA)-0.4295*log10(T_c_steady(i,1,AA)));                                                    %mu matrix   
end


%% Initial Force condition acting on the FAST
% Buoyancy
for i=1:node_MATRA-1
        alpha=floor(zi_transient_F(1,i,AA)/dz)+1;                   %coolant node 의 정수
        gamma=zi_transient_F(1,i,AA)/dz;
        beta=floor(zi_transient_F(1,i,AA)/dz);                      %coolant node 의 정수
        if beta <=0
            beta=beta+1;
            alpha=alpha+1;
            gamma=gamma+1;
        end
        if i==1
            Fb(i,AA)=-(den_c(beta)*(alpha-gamma)+den_c(alpha)*(gamma-beta))*g*A_front*dh*0.5;
        end
        if i==node_MATRA-1
            Fb(i,AA)=-(den_c(beta)*(alpha-gamma)+den_c(alpha)*(gamma-beta))*g*A_front*dh*0.5;
        end
        Fb(i,AA)=-(den_c(beta)*(alpha-gamma)+den_c(alpha)*(gamma-beta))*g*A_front*dh;
end
F_buoyancy(AA)=sum(Fb(:,AA));
%gravity
F_gravity =density_FAST* V * g ;
%drag and pressure
F_drag(AA)=0;     F_pressure(AA)=0;
%% Basic properties
Q_error(AA)=0;
Pe_dz(AA)=0;
a(AA)=0;
v(AA)=a(AA)*dt;
s(AA)=v(AA)*dt;                              
inserted_reactivity_FAST = 0;        %initial inserted reactivity
feedback_reactivity      = 0;        %initial feedback reactivity


%% triger condition ( the condition that FAST could fall down)
trigger(AA)=F_gravity+F_buoyancy(AA);


dt=0.01;
control_rod_withrawal=68.7*cent;
% store steady-state fuel and coolant temperature
Tf_mean_steady=mean(mean(mean(Plot_TC(1:node_f,1:end,:))));                        %fuel mean temperature at steady-state
Tc_mean_steady=mean(mean(Plot_TC(node_t+1,1:end,:)));                        %coolatn mean temperature at steady-state
Tf_mean_transient=mean(mean(Plot_TC(1:node_f,1:end,:)));                     %fuel mean temperature at transient
Tc_mean_transient=mean(mean(Plot_TC(node_t+1,1:end,:)));                     %coolatn mean temperature at transient

transient_condition=1;                                                     %출력 증가 상황
iteration=1;        %5000마다 리셋
iteration_main=1;
iteration_cut=2000; %5000마다 데이터 저장
er=1;
tol2=1e-6;
q_t(AA)=q*transient_condition*Neutron(iteration)*Group_Q(AA);                                                                   %heat generation 392.
v_i(:,:,AA)=zeros(num_gap,1);                                                      %velocity distribution
v_i_p(:,:,AA)=zeros(num_gap,1);         
con_q(AA)=q_t(AA)/(pi*f_o^2*(2*H_a/pi));                                                %0.63662는 cos를 높이로 적분한 값에 비롯됨. 77page 연구노트

    %previous velocity distribution
end
% pressure gradient algorithm
a_coef=0.0001;
Pe_dz_a=30;
Pe_dz_b=-30;
tol2=1e-6; % for pressure gradient
Pe_dz=zeros(Group_N,1);
Time_limit=80;
tic
while Ldt<=Time_limit                     %Trainsient 때 사용함

for AA = 1: length(Group_Q)
%power

con_q(AA)=q_t(AA)/(pi*f_o^2*(2*H_a/pi));                                                %0.63662는 cos를 높이로 적분한 값에 비롯됨. 77page 연구노트

for j=1:node_t                                                              %for overall region including fuel and clad
    for i=1:node_MATRA
        % volume & area
        if j < node_f                                                      % Fuel 영역에서
            V_rf(i,j)=2*pi*ddfh*ddf*n_t(j);                                % unit volume
            A_rf(i,j)=2*pi*(n_t(j)+ddf/2)*ddfh;                            % node right side area
            A_lf(i,j)=2*pi*(n_t(j)-ddf/2)*ddfh;                            % node left side area
            A_uf(i,j)=2*pi*ddf*n_t(j);                                     % node up side area
            A_df(i,j)=2*pi*ddf*n_t(j);                                     % node up side area
            if j==1                                                        % 왼쪽 단열 (중심 부에서)
                V_rf(i,j)=pi/4*ddf^2*ddfh;
                A_lf(i,j)=0;                         
                A_uf(i,j)=pi/4*ddf^2;
                A_df(i,j)=A_uf(i,j);
            end

            if i==1                                                             %맨 윗면에 대해서 면적 및 부피가 기존것의 반. dh가 절반이기 때문
                V_rf(i,j)=V_rf(i,j)/2;
                A_rf(i,j)=A_rf(i,j)/2;
                A_lf(i,j)=A_lf(i,j)/2;
                A_uf(i,j)=0;                                                     %위쪽방향으로 열전달 x 
            end
            if i==node_MATRA                                                    %맨 아랫면에 대해서 면적 및 부피가 기존의 반. dh가 절반이기 때문
                V_rf(i,j)=V_rf(i,j)/2;
                A_rf(i,j)=A_rf(i,j)/2;
                A_lf(i,j)=A_lf(i,j)/2;
                A_df(i,j)=0;                                                     %아래방향으로 열전달 x
            end 
        end


        if j > node_f                                                      %Clad 영역에서
            V_rf(i,j)=2*pi*ddfh*ddc*n_t(j);                                %unit volume
            A_rf(i,j)=2*pi*(n_t(j)+ddc/2)*ddfh;                            %node right side area
            A_lf(i,j)=2*pi*(n_t(j)-ddc/2)*ddfh;                            %node left side area
            A_uf(i,j)=2*pi*ddc*n_t(j);                                     %node up side area
            A_df(i,j)=2*pi*ddc*n_t(j);                                     %node up side area

            if j==node_t                                                   %맨 오른쪽 경계면에서 ( 반쪼가리)
                V_rf(i,j)=pi*ddfh*ddc*(n_t(j)-ddc/4);
                A_rf(i,j)=2*pi*(n_t(j))*ddfh;
                A_uf(i,j)=pi*ddc*(n_t(j)-ddc/4);
                A_df(i,j)=A_uf(i,j);
            end
            if i==1                                                        %맨 윗면에 대해서 면적 및 부피가 기존것의 반. dh가 절반이기 때문
                V_rf(i,j)=V_rf(i,j)/2;
                A_rf(i,j)=A_rf(i,j)/2;
                A_lf(i,j)=A_lf(i,j)/2;
                A_uf(i,j)=0;                                               %위쪽방향으로 열전달 x 
            end
            if i==node_MATRA                                               %맨 아랫면에 대해서 면적 및 부피가 기존의 반. dh가 절반이기 때문
                V_rf(i,j)=V_rf(i,j)/2;
                A_rf(i,j)=A_rf(i,j)/2;
                A_lf(i,j)=A_lf(i,j)/2;
                A_df(i,j)=0;                                               %아래방향으로 열전달 x
            end 
        end


        if j == node_f                                                         %경계면에서 Fuel and Cladding interface
            V_ff = pi*ddfh*ddf*(n_t(j)-ddf/4) ;                               %volume fraction of fuel
            V_fc = pi*ddfh*ddc*(n_t(j)+ddc/4) ;                                %volume fraction of cla
            
            V_rf(i,j)=V_ff+V_fc;    %unit volume
            A_rf(i,j)=2*pi*(n_t(j)+ddc/2)*ddfh;                                 %node right side area
            A_lf(i,j)=2*pi*(n_t(j)-ddf/2)*ddfh;                                 %node left side area
            A_uf(i,j)=pi*ddf*(n_t(j)-ddf/4)+pi*ddc*(n_t(j)+ddc/4);              %node up side area
            A_df(i,j)=pi*ddf*(n_t(j)-ddf/4)+pi*ddc*(n_t(j)+ddc/4);              %node up side area
            if i==1                                                             %맨 윗면에 대해서 면적 및 부피가 기존것의 반. dh가 절반이기 때문
                V_rf(i,j)=V_rf(i,j)/2;
                A_rf(i,j)=A_rf(i,j)/2;
                A_lf(i,j)=A_lf(i,j)/2;
                A_uf(i,j)=0;                                                     %위쪽방향으로 열전달 x 
            end
            if i==node_MATRA                                                    %맨 아랫면에 대해서 면적 및 부피가 기존의 반. dh가 절반이기 때문
                V_rf(i,j)=V_rf(i,j)/2;
                A_rf(i,j)=A_rf(i,j)/2;
                A_lf(i,j)=A_lf(i,j)/2;
                A_df(i,j)=0;                                                     %아래방향으로 열전달 x
            end 
        end        

        %%properties
        %conductivity, density, capacity
        if j<=node_f                                      %Uzr region
            k_uzr(i,j)=S(T_uzr((i-1)*node_t+j,1,AA));     %[W/m.K]
            den_uzr(i,j)=12375.864; %density matrix
            c_uzr(i,j)=121.336;                                        %capacity matrix
        else                                            %HT-9 region
            k_uzr(i,j)=17.622+2.428*10^(-2)*(T_uzr((i-1)*node_t+j,1,AA))-1.696*10^(-5)*(T_uzr((i-1)*node_t+j,1,AA))^2;%[W/m.K]
            if  T_uzr((i-1)*node_t+j,1,AA)<800
                c_uzr(i,j)=561;%(T_uzr((i-1)*node_t+j,1)-500)/6+500;
            else
                c_uzr(i,j)=561;%(T_uzr((i-1)*node_t+j,1)-800)/6+550;
            end
             den_uzr(i,j)=7750;% 7700~7800 ferritic S/S 
        end
       
        %properties at the interface follows volume fraction.
        if j==node_f
            k_uzr(i,j)=(S(T_uzr((i-1)*node_t+j,1,AA)) + 17.622+2.428*10^(-2)*(T_uzr((i-1)*node_t+j,1,AA))-1.696*10^(-5)*(T_uzr((i-1)*node_t+j,1,AA))^2 )/2;      %[W/m.K]

        end
        
        
    end
end

%% 좌변/우변에 곱해지는 특성치들을 고려하기 위해서   KL = 왼쪽 노드, KR= 오른쪽 노드
for i=1:node_MATRA
    for j=1:node_t
        if j<=node_f
            KL(i,j)=S(T_uzr((i-1)*node_t+j,1,AA));                                %[W/m.K] Fuel side
        else
            KL(i,j)=17.622+2.428*10^(-2)*(T_uzr((i-1)*node_t+j,1,AA))-1.696*10^(-5)*(T_uzr((i-1)*node_t+j,1,AA))^2;   %[W/m.K] cladding side
        end
    end
    for j=1:node_t
        if j<=node_f-1
            KR(i,j)=S(T_uzr((i-1)*node_t+j,1,AA));                                %[W/m.K] Fuel side
        else
            KR(i,j)=17.622+2.428*10^(-2)*(T_uzr((i-1)*node_t+j,1,AA))-1.696*10^(-5)*(T_uzr((i-1)*node_t+j,1,AA))^2;   %[W/m.K] cladding side
        end    
    end
end

%% 조화 평균 conductivity
% KL_m=KL; % 변수 선언
% KR_m=KR; % 변수 선언
% % Centerline에 대해서
% KL_m(1:end,1)=0;                                                           % 중앙에서 왼쪽으로 열전달 x
% KR_m(1:end,1)=k_uzr(1:end,1);                                              % 오른쪽으로 열 전달시 conductivity는 centerline것을 이용. ( 추후 변경 가능 k i+1/2 로)
% % 조화평균 사용                                                             %
% % km(2*k1k2/(k1+k2)
% for i=2:length(KL(1,1:end))
%     k1=KL(1:end,i-1);
%     k2=KL(1:end,i);
%     KL_m(1:end,i)=2*k1.*k2./(k1+k2);                                       % 좌변 노드에 대한 물성치
%     
%     k1=KR(1:end,i-1);
%     k2=KR(1:end,i);
%     KR_m(1:end,i)=2*k1.*k2./(k1+k2);                                       % 우변 노드에 대한 물성치
% end
%Conductivity
k_uzr_m=zeros(length(k_uzr(1:end,1)),length(k_uzr(1,1:end))-1);            %매트릭스 생성
for i=2:length(KL(1,1:end))
    k1=k_uzr(1:end,i-1);
    k2=k_uzr(1:end,i);
    k_uzr_m(1:end,i-1)=2*k1.*k2./(k1+k2); 
end
%중심에서 열전도계수
k_uzr_m=[k_uzr(1:end,1),k_uzr_m];

%heat flux
A_q=A_df(2,1:node_f);  A_q(1,node_f)=pi*ddf*(n_t(node_f)-ddf/4);                           %heat flux에 곱해질 수
n_h=[H_a/2,(H_a/2-ddfh/2):-ddfh:(-H_a/2+ddfh/2),-H_a/2];                    %heat flux에 적분될 변수들
% q_m = heat generation
for i=1:node_MATRA
    for j=1:node_f                                                         %for overall region including fuel and clad
        q_m(i,j)=A_q(j)/pi*con_q(AA)*H_a*(sin(pi/H_a*n_h(i))-sin(pi/H_a*n_h(i+1)))  ;  % heat generation node
    end
end
% q_mv=volumetric heat generation
q_mv=q_m./V_rf;


%%Active core region coolant side
for i=1:node_MATRA
%% liquid sodium properties
    k_coolant(i,1)=92.951-5.8087e-2*(T_coolant(i,1,AA)-273.15)+11.7274e-6*(T_coolant(i,1,AA)-273.15)^2;                                       %k matrix
    den_coolant(i,1,AA)=950.0483-0.2298*(T_coolant(i,1,AA)-273.15)-14.6045e-6*(T_coolant(i,1,AA)-273.15)^2+5.6377e-9*(T_coolant(i,1,AA)-273.15)^3; %density matrix
    cp_coolant(i,1)=  1436.715-0.5806*(T_coolant(i,1,AA)-273.15)+4.6273e-4*(T_coolant(i,1,AA)-273.15)^2; %          Coolant_cp(T_coolant(i,1)-273.15); %                           %capacity matrix
    mu_coolant(i,1)=10^(-2.4892+220.65/T_coolant(i,1,AA)-0.4295*log10(T_coolant(i,1,AA)));                                                    %mu matrix   
    Re(i,1,AA)=G*D_fh/mu_coolant(i,1);                                                                                              %Reynolds number
    Pr(i,1)=cp_coolant(i,1)*mu_coolant(i,1)/k_coolant(i,1);                                                                                         %Plantle number
    H_c(i,1)=k_coolant(i,1)/D_fh*0.023*Re(i,1)^0.8*Pr(i,1)^0.4;                                                                   %Dittus-Boelter equation. heat transfer coefficient

end
for j=1:node_r
    for i=1:num_h
        % volume & area
        V_r(i,j)=2*pi*ddh*ddr*n_r(j);           %unit volume
        A_r(i,j)=2*pi*(n_r(j)+ddr/2)*ddh;       %node right side area
        A_l(i,j)=2*pi*(n_r(j)-ddr/2)*ddh;       %node left side area
        A_u(i,j)=2*pi*ddr*n_r(j);               %node up side area
        A_d(i,j)=2*pi*ddr*n_r(j);               %node up side area
        if j==1
            V_r(i,j)=pi/4*ddr^2*ddh;
            A_l(i,j)=0;                         %왼쪽 단열
            A_u(i,j)=pi/4*ddr^2;
            A_d(i,j)=A_u(i,j);
        end
        if j==node_r
            V_r(i,j)=pi*ddh*ddr*(n_r(j)-ddr/4);
            A_r(i,j)=2*pi*(n_r(j))*ddh;
            A_u(i,j)=pi*ddr*(n_r(j)-ddr/4);
            A_d(i,j)=A_u(i,j);
        end
        if i==1
            V_r(i,j)=V_r(i,j)/2;
            A_r(i,j)=A_r(i,j)/2;
            A_l(i,j)=A_l(i,j)/2;
            A_u(i,j)=0;
        end
        if i==num_h
            V_r(i,j)=V_r(i,j)/2;
            A_r(i,j)=A_r(i,j)/2;
            A_l(i,j)=A_l(i,j)/2;
            A_d(i,j)=0;
        end 
         % properties
         %conductivity
        if j<=floor((D_p/2-t_cf)/ddr)                    %fluid region
            k_f(i,j)=92.951-5.8087e-2*(T_f((i-1)*node_r+j,1,AA)-273.15)+11.7274e-6*(T_f((i-1)*node_r+j,1,AA)-273.15)^2;     %[W/m.K]
            den_f(i,j)=950.0483-0.2298*(T_f((i-1)*node_r+j,1,AA)-273.15)-14.6045e-6*(T_f((i-1)*node_r+j,1,AA)-273.15)^2+5.6377e-9*(T_f((i-1)*node_r+j,1,AA)-273.15)^3; %density matrix
            c_f(i,j)=1436.715-0.5806*(T_f((i-1)*node_r+j,1,AA)-273.15)+4.6273e-4*(T_f((i-1)*node_r+j,1,AA)-273.15)^2;    %Coolant_cp((T_f((i-1)*node_r+j,1)-273.15));% 
            
            %                                   %capacity matrix
        else                                            %HT-9 region
            k_f(i,j)=17.622+2.428*10^(-2)*(T_f((i-1)*node_r+j,1,AA))-1.696*10^(-5)*(T_f((i-1)*node_r+j,1,AA))^2;%[W/m.K]
            if T_f((i-1)*node_r+j,1)<800
                c_f(i,j)=(T_f((i-1)*node_r+j,1,AA)-500)/6+500;
            else
                c_f(i,j)=(T_f((i-1)*node_r+j,1,AA)-800)/6+550;
            end
             den_f(i,j)=7750;% 7700~7800 ferritic S/S 
        end
    end
end

%% Fuel and cladding Matrix
mcf =den_uzr.*V_rf.*c_uzr./dt;                                             %Matrix coefficient (alpha)
% boundary condition at the interface
density_fuelside= 12375.864*ones(node_MATRA,1);     %[kg/m3]
density_cladside= 7750*ones(node_MATRA,1);% 
volume_fuelside = pi*ddfh*ddf*(n_t(node_f)-ddf/4)*ones(node_MATRA,1) ;                               %volume fraction of fuel
volume_cladside = pi*ddfh*ddc*(n_t(node_f)+ddc/4)*ones(node_MATRA,1) ;                               %volume fraction of cla
volume_fuelside(1,1)=volume_fuelside(1,1)/2; volume_fuelside(end,1)=volume_fuelside(end,1)/2;        %end point is half
volume_cladside(1,1)=volume_cladside(1,1)/2; volume_cladside(end,1)=volume_cladside(end,1)/2;
capacity_fuelside=121.336;
capacity_cladside=561;

%pcp/dt*V
mcf(1:end,node_f)=(density_fuelside.*volume_fuelside*capacity_fuelside./dt)+(density_cladside.*volume_cladside*capacity_cladside./dt);
%pcp/dt
vcf=mcf./V_rf;

%% FAST region matrix
mc =den_f.*V_r.*c_f./dt;            
vc =mc./V_r;


%% matrix 생성
% Interior node 에 대해서 먼저 실행
NN=node_t;
MM=node_MATRA;                                % z 방향 노드
FF=node_r;                                    % FAST r 방향 노드
CN=1;                                         % 쿨런트 r 방향 노드
TT=NN+FF+CN;                                  % 총 노드




%%  Interior node
for i=1*TT+1:(MM-1)*TT-1
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end

    % FUEL ROD REGION
    if (node_j>1) && (node_j<NN)
        %heat generation
        MH(i,1)=q_mv(node_i,node_j)/vcf(node_i,node_j);
        % Tob ( conductivity 는 중앙 노드꺼씀)
        TC=-k_uzr(node_i,node_j)/ddfh^2/vcf(node_i,node_j);                        
        % BOT ( conductivity 는 중앙 노드꺼씀)
        BC=-k_uzr(node_i,node_j)/ddfh^2/vcf(node_i,node_j);                         
        % Left (조화평균값)
        LC=-n_frh(node_j)*k_uzr_m(node_i,node_j)    /n_t(node_j)/(n_frh(node_j+1)-n_frh(node_j))/vcf(node_i,node_j)/(n_t(node_j)  -n_t(node_j-1));
        % Right (conductivity 조화 평균)
        RC=-n_frh(node_j+1)*k_uzr_m(node_i,node_j+1)/n_t(node_j)/(n_frh(node_j+1)-n_frh(node_j))/vcf(node_i,node_j)/(n_t(node_j+1)-n_t(node_j));
        % Center
        CC=-(TC+BC+LC+RC)+1;

        M(i,i-TT)=TC;
        M(i,i+1)=RC;
        M(i,i)=CC;
        M(i,i-1)=LC;
        M(i,i+TT)=BC;
    
    end
    
        % FAST region
    if (node_j>NN+CN+1) && (node_j<TT)
        node_j=node_j-(NN+CN);                                             %인덱스 조정
        % Tob ( conductivity 는 중앙 노드꺼씀)
        TC=-k_f(node_i,node_j)/ddh^2/vc(node_i,node_j);                        
        % BOT ( conductivity 는 중앙 노드꺼씀)
        BC=-k_f(node_i,node_j)/ddh^2/vc(node_i,node_j);                         
        % Left (조화평균값)
        LC=-n_rF(node_j)*k_f(node_i,node_j)    /n_r(node_j)/(n_rF(node_j+1)-n_rF(node_j))/vc(node_i,node_j)/(n_r(node_j)  -n_r(node_j-1));
        % Right (conductivity 조화 평균)
        RC=-n_rF(node_j+1)*k_f(node_i,node_j)/n_r(node_j)/(n_rF(node_j+1)-n_rF(node_j))/vc(node_i,node_j)/(n_r(node_j+1)-n_r(node_j));
        % Center
        CC=-(TC+BC+LC+RC)+1;

        M(i,i-TT)=TC;
        M(i,i+1)=RC;
        M(i,i)=CC;
        M(i,i-1)=LC;
        M(i,i+TT)=BC;
    
    end
    
    
end
%%  Fuel centerline node
for i=1:TT:TT*MM
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    %heat generation
    MH(i,1)=q_mv(node_i,node_j)/vcf(node_i,node_j);

    % Right (conductivity 조화 평균)
    RC=-4*k_uzr_m(node_i,node_j)/vcf(node_i,node_j)/(n_t(node_j+1)-n_t(node_j))^2;
    % Center
    CC=-(RC)+1;
    
    M(i,i+1)=RC;
    M(i,i)=CC;
end
%% FAST centerline node
for i=NN+CN+1:TT:TT*MM
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    node_j=node_j-(NN+CN);                                                 %인덱스 조정
    %  (conductivity 일반 값 사용)
    RC=-4*k_f(node_i,node_j)/vc(node_i,node_j)/(n_r(node_j+1)-n_r(node_j))^2;
    % Center
    CC=-(RC)+1;
    
    M(i,i+1)=RC;
    M(i,i)=CC;
    
end

%%  Fuel upper node
for i=2:NN-1
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    %heat generation
    MH(i,1)=q_mv(node_i,node_j)/vcf(node_i,node_j);

    % BOT ( conductivity 는 중앙 노드꺼씀)
    BC=-k_uzr(node_i,node_j)/ddfh^2/vcf(node_i,node_j);                         
    % Left (조화평균값)
    LC=-n_frh(node_j)*k_uzr_m(node_i,node_j)/n_t(node_j)/(n_frh(node_j+1)-n_frh(node_j))/vcf(node_i,node_j)/(n_t(node_j)-n_t(node_j-1));
    % Right (conductivity 조화 평균)
    RC=-n_frh(node_j+1)*k_uzr_m(node_i,node_j+1)/n_t(node_j)/(n_frh(node_j+1)-n_frh(node_j))/vcf(node_i,node_j)/(n_t(node_j+1)-n_t(node_j));
    % Center
    CC=-(BC+LC+RC)+1;
  

    M(i,i+1)=RC;
    M(i,i)=CC;
    M(i,i-1)=LC;
    M(i,i+TT)=BC;
    
end
%%  FAST upper node
for i=NN+CN+2:TT-1
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    node_j=node_j-(NN+CN);                                             %인덱스 조정
    % BOT ( conductivity 는 중앙 노드꺼씀)
    BC=-k_f(node_i,node_j)/ddh^2/vc(node_i,node_j);                         
    % Left (조화평균값)
    LC=-n_rF(node_j)*k_f(node_i,node_j)/n_r(node_j)/(n_rF(node_j+1)-n_rF(node_j))/vc(node_i,node_j)/(n_r(node_j)-n_r(node_j-1));
    % Right (conductivity 조화 평균)
    RC=-n_rF(node_j+1)*k_f(node_i,node_j)/n_r(node_j)/(n_rF(node_j+1)-n_rF(node_j))/vc(node_i,node_j)/(n_r(node_j+1)-n_r(node_j));
    % Center
    CC=-(BC+LC+RC)+1;
  

    M(i,i+1)=RC;
    M(i,i)=CC;
    M(i,i-1)=LC;
    M(i,i+TT)=BC;
    
end

%%  Fuel bottom node
for i=(MM-1)*(TT)+2:MM*TT-(FF+2)
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    %heat generation
    MH(i,1)=q_mv(node_i,node_j)/vcf(node_i,node_j);
    % Tob ( conductivity 는 중앙 노드꺼씀)
    TC=-k_uzr(node_i,node_j)/ddfh^2/vcf(node_i,node_j);   
    % Left (조화평균값)
    LC=-n_frh(node_j)*k_uzr_m(node_i,node_j)    /n_t(node_j)/(n_frh(node_j+1)-n_frh(node_j))/vcf(node_i,node_j)/(n_t(node_j)  -n_t(node_j-1));
    % Right (conductivity 조화 평균)
    RC=-n_frh(node_j+1)*k_uzr_m(node_i,node_j+1)/n_t(node_j)/(n_frh(node_j+1)-n_frh(node_j))/vcf(node_i,node_j)/(n_t(node_j+1)-n_t(node_j));
    % Center
    CC=-(TC+LC+RC)+1;
  
    M(i,i-TT)=TC;
    M(i,i+1)=RC;
    M(i,i)=CC;
    M(i,i-1)=LC;
end
%%  FAST bottom node
for i=(MM-1)*(TT)+(NN+CN)+2:MM*TT-(1)
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    node_j=node_j-(NN+CN);                                             %인덱스 조정

    % Tob ( conductivity 는 중앙 노드꺼씀)
    TC=-k_f(node_i,node_j)/ddh^2/vc(node_i,node_j);   
    % Left (조화평균값)
    LC=-n_rF(node_j)*k_f(node_i,node_j)    /n_r(node_j)/(n_rF(node_j+1)-n_rF(node_j))/vc(node_i,node_j)/(n_r(node_j)  -n_r(node_j-1));
    % Right (conductivity 조화 평균)
    RC=-n_rF(node_j+1)*k_f(node_i,node_j)/n_r(node_j)/(n_rF(node_j+1)-n_rF(node_j))/vc(node_i,node_j)/(n_r(node_j+1)-n_r(node_j));
    % Center
    CC=-(TC+LC+RC)+1;
  
    M(i,i-TT)=TC;
    M(i,i+1)=RC;
    M(i,i)=CC;
    M(i,i-1)=LC;
end

%%  Fuel rod outter node
for i=NN:TT:TT*MM
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end
    %heat generation
    MH(i,1)=q_mv(node_i,node_j)/vcf(node_i,node_j);

    % Tob ( conductivity 는 중앙 노드꺼씀)
    TC=-k_uzr(node_i,node_j)/ddfh^2/vcf(node_i,node_j);                        
    % BOT ( conductivity 는 중앙 노드꺼씀)
    BC=-k_uzr(node_i,node_j)/ddfh^2/vcf(node_i,node_j);                         
    % Left (조화평균값)
    LC=-k_uzr_m(node_i,node_j)*A_lf(node_i,node_j)/(n_t(node_j)-n_t(node_j-1))/mcf(node_i,node_j);
    % Right (conductivity 조화 평균)
    RC=-H_c(node_i,1)*A_rf(node_i,node_j)/mcf(node_i,node_j);
    CC=-(TC+BC+LC+RC)+1;

        if node_i==1
            CC=-(BC+LC+RC)+1;                                                                   % 맨 위쪽 열전달 x
            M(i,i+1)=RC;
            M(i,i)=CC;
            M(i,i-1)=LC;
            M(i,i+TT)=BC;
        end
        if (1<node_i) && (node_i<MM)
            CC=-(TC+BC+LC+RC)+1;
            M(i,i-TT)=TC;
            M(i,i+1)=RC;
            M(i,i)=CC;
            M(i,i-1)=LC;
            M(i,i+TT)=BC;
        end
        if node_i==MM                                                        % 맨 아래쪽 열전달 x
            CC=-(TC+LC+RC)+1;
            M(i,i-TT)=TC;
            M(i,i+1)=RC;
            M(i,i)=CC;
            M(i,i-1)=LC;
        end 

    

end

%%  FAST pin outter node
for i=TT:TT:TT*MM
    node_j=rem(i,TT);                  % node j index
    node_i=fix(i/TT)+1;                % node i index
    if node_j==0
        node_j=TT;
        node_i=node_i-1;
    end

    node_j=node_j-(NN+CN);                                             %인덱스 조정

    % Tob ( conductivity 는 중앙 노드꺼씀)
    TC=-k_f(node_i,node_j)/ddh^2/vc(node_i,node_j);                        
    % BOT ( conductivity 는 중앙 노드꺼씀)
    BC=-k_f(node_i,node_j)/ddh^2/vc(node_i,node_j);                         
    % Left (조화평균값)
    LC=-k_f(node_i,node_j)*A_l(node_i,node_j)/(n_r(node_j)-n_r(node_j-1))/mc(node_i,node_j);
    % Right (conductivity 조화 평균)
    RC=-H_c(node_i,1)*A_r(node_i,node_j)/mc(node_i,node_j);
    CC=-(TC+BC+LC+RC)+1;

        if node_i==1
            CC=-(BC+LC+RC)+1;                                                                   % 맨 위쪽 열전달 x
            M(i,i-FF)=RC;
            M(i,i)=CC;
            M(i,i-1)=LC;
            M(i,i+TT)=BC;
        end
        if (1<node_i) && (node_i<MM)
            CC=-(TC+BC+LC+RC)+1;
            M(i,i-TT)=TC;
            M(i,i-FF)=RC;
            M(i,i)=CC;
            M(i,i-1)=LC;
            M(i,i+TT)=BC;
        end
        if node_i==MM                                                        % 맨 아래쪽 열전달 x
            CC=-(TC+LC+RC)+1;
            M(i,i-TT)=TC;
            M(i,i-FF)=RC;
            M(i,i)=CC;
            M(i,i-1)=LC;
        end 

end

%% Coolant hyperbolic equation
for nodei=NN+CN:TT:TT*MM
    j=rem(nodei,TT);                  % node j index
    i=fix(nodei/TT)+1;                % node i index
    if j==0
        j=TT;
        i=i-1;
    end
%% Hyperboilic equation (Coolant matrix)
a1=Massflux_FA(AA)*dt/den_coolant(i,1,AA)/(n_h(1,i)-n_h(1,i+1));                            %Coolant in and out energy transfer
b1=1/2*A_rf(i,end)*dt/den_coolant(i,1,AA)/(A_c_f)/(n_h(1,i)-n_h(1,i+1))/cp_coolant(i,1)*H_c(i,1);  %heat source-coolant interaction ( 1/3 = 120/360)
c1=1/6*A_rf(i,end)*dt/den_coolant(i,1,AA)/(A_c_f)/(n_h(1,i)-n_h(1,i+1))/cp_coolant(i,1)*H_c(i,1);  %fast-coolant interaction (1/6 = 60/360 )
M(nodei,nodei)=1+a1+b1+c1;                                                 % - center matrix variable
if i<MM
M(nodei,nodei+TT)=-a1;                                                     % - previous coolant matrix variable
else
MH(nodei,1)=a1*T_inlet;                                                   %inlet boundary
end
M(nodei,nodei-1)=-b1;                                                      % - cladding - coolant
M(nodei,nodei+FF)=-c1;                                                     % - fast - coolant

end

T_next(:,:,AA)=M\T_T(:,:,AA)+M\MH;




for i=1:node_MATRA
    for j=1:node_total
        if j<=node_t
            T_uzr(j+(i-1)*node_t,1,AA)=T_T(  node_total*(i-1)+j ,1,AA);
        end
        if j==node_t+1
            T_coolant(i,1,AA)=T_T(  node_total*(i-1)+j ,1,AA);
    
        end
        if j>node_t+1
            T_f((j-(node_t+1))+(i-1)*node_r,1,AA)=T_T(  node_total*(i-1)+j ,1,AA);
        end
    end
end
T_T(:,:,AA)=T_next(:,:,AA);
Plot_T(:,:,AA)=reshape(T_next(:,:,AA),node_total,node_MATRA);
Plot_TA(:,:,AA)=Plot_T(1:node_t+1,:,AA);
Plot_TB(:,:,AA)=flipud(Plot_T(node_t+2:end,:,AA));
Plot_TC(:,:,AA)=[Plot_TA(:,:,AA);Plot_TB(:,:,AA)];
%    surf(Plot_TC)
%  AAA=sum((Plot_TC(node_t,1:end)-Plot_TC(node_t+1,1:end)).*(reshape(H_c,1,node_MATRA)).*reshape(A_r(1:end,end),1,node_MATRA))
%  BBB=G*A_c_f*2*cp_coolant(1,1)*(Plot_TC(node_t+1,1)-Plot_TC(node_t+1,end))




Interpx=H_a:-ddfh:0;                                %인터폴레이션 x 축
Interpy=Plot_TC(end-round(R/ddr),1:end,AA);            %FAST 벽면에 해당하는 온도    
for i=1:node_MATRA                              %interpolation을 통해 온도 매칭.
    T_c_steady(i,1,AA)=interp1(Interpx,Interpy,(i-1)*dz);
    if (i-1)*dz>H_a
        T_c_steady(i,1,AA)=T_c_steady(i-1,1,AA);
    end       
end

for i=1:node_MATRA
%liquid sodium properties
    k_c(i,1)=92.951-5.8087e-2*(T_c_steady(i,1,AA)-273.15)+11.7274e-6*(T_c_steady(i,1,AA)-273.15)^2;                                       %k matrix
    den_c(i,1)=950.0483-0.2298*(T_c_steady(i,1,AA)-273.15)-14.6045e-6*(T_c_steady(i,1,AA)-273.15)^2+5.6377e-9*(T_c_steady(i,1,AA)-273.15)^3; %density matrix
    cp_c(i,1)=1436.715-0.5806*(T_c_steady(i,1,AA)-273.15)+4.6273e-4*(T_c_steady(i,1,AA)-273.15)^2;                                        %capacity matrix
    mu_c(i,1)=10^(-2.4892+220.65/T_c_steady(i,1,AA)-0.4295*log10(T_c_steady(i,1,AA)));                                                    %mu matrix   
end



if v(AA)>0
    reverse=-1;
else
    reverse=1;
end

%% 위로 올라가려는 조건
if zi_transient_F(1,end,AA) >= H_p                                          %벽에 다시 닿고 힘을 위로 받을 때
    %% each time step data saving in steady state
    a(AA)=0;
    v(AA)=0;
    s(AA)=0;
    F_drag(AA)=0;
    F_pressure(AA)=0;
    Q_error(AA)=0;
    Pe_dz(AA)=0;
    zi_transient_F(:,:,AA)=zi_F+(H_p-H);

        
    plota(iteration,1,AA)=a(AA);
    plotv(iteration,1,AA)=v(AA);
    plots(iteration,1,AA)=s(AA);
    plotb(iteration,1,AA)=F_buoyancy(AA);
    plotg(iteration,1,AA)=F_gravity;
    plott(iteration,1,AA)=trigger(AA);
    plotd(iteration,1,AA)=F_drag(AA);
    plotp(iteration,1,AA)=F_pressure(AA);
    plotq(iteration,1,AA)=abs(v(AA))*A_front;                             %mass flow FAST pusing
    plotpd(iteration,1,AA)=Pe_dz(AA);
      %  zi_transient=zi_transient+s;
      %actual FAST position and each meshing point
    plot_reactivity(iteration,1)=total_reactivity;
    plot_feedback(iteration,1)=feedback_reactivity;
    plot_inserted_r(iteration,1)=inserted_reactivity_FAST;
    plot_position_start(iteration,1,AA)=zi_transient_F(1,1,AA);
    plot_position_end(iteration,1,AA)=zi_transient_F(1,end,AA);
    plot_absorber_h(iteration,1,AA)=ab_transient_F+plots(iteration,1,AA);                         % time-dependent absorber position
end



%% bouyancy
% Buoyancy
for i=1:node_MATRA-1
        alpha=floor(zi_transient_F(1,i,AA)/dz)+1;                   %coolant node 의 정수
        gamma=zi_transient_F(1,i,AA)/dz;
        beta=floor(zi_transient_F(1,i,AA)/dz);                      %coolant node 의 정수
        if beta<1 
            alpha=alpha+1;
            gamma=gamma+1;
            beta=beta+1;
        end
        if i==1
            Fb(i,AA)=-(den_c(beta)*(alpha-gamma)+den_c(alpha)*(gamma-beta))*g*A_front*dh*0.5;
            plot_den_b(iteration,AA)=(den_c(beta)*(alpha-gamma)+den_c(alpha)*(gamma-beta));

        end
        if i==node_MATRA-1
            Fb(i,AA)=-(den_c(beta)*(alpha-gamma)+den_c(alpha)*(gamma-beta))*g*A_front*dh*0.5;
            plot_den_t(iteration,AA)=(den_c(beta)*(alpha-gamma)+den_c(alpha)*(gamma-beta));

        end
        Fb(i,AA)=-(den_c(beta)*(alpha-gamma)+den_c(alpha)*(gamma-beta))*g*A_front*dh;
end
xe = eq_point(dz,node_MATRA,dh,g,A_front,den_c,density_FAST,V,H,H_p);
plot_xe(iteration,Group_N)=xe; 

F_buoyancy(AA)=sum(Fb(:,AA));
F_gravity =density_FAST* V * g ;
%% triger condition ( the condition that FAST could fall down)
trigger(AA)=F_gravity+F_buoyancy(AA);

%% volumetric flow calculation  < inlet condition has free developed velocity profile >
if trigger(AA)>=0 && zi_transient_F(1,end,AA)==H_p
        %% each time step data saving in steady state
    a(AA)=0;
    v(AA)=0;
    s(AA)=0;
    F_drag(AA)=0;
    F_pressure(AA)=0;
    Q_error(AA)=0;
    Pe_dz(AA)=0;
    zi_transient_F(:,:,AA)=zi_F+(H_p-H);

        
    plota(iteration,1,AA)=a(AA);
    plotv(iteration,1,AA)=v(AA);
    plots(iteration,1,AA)=s(AA);
    plotb(iteration,1,AA)=F_buoyancy(AA);
    plotg(iteration,1,AA)=F_gravity;
    plott(iteration,1,AA)=trigger(AA);
    plotd(iteration,1,AA)=F_drag(AA);
    plotp(iteration,1,AA)=F_pressure(AA);
    plotq(iteration,1,AA)=abs(v(AA))*A_front;                             %mass flow FAST pusing
      %  zi_transient=zi_transient+s;
      %actual FAST position and each meshing point
    plot_reactivity(iteration,1)=total_reactivity;
    plot_feedback(iteration,1)=feedback_reactivity;
    plot_inserted_r(iteration,1)=inserted_reactivity_FAST;
    plot_position_start(iteration,1,AA)=zi_transient_F(1,1,AA);
    plot_position_end(iteration,1,AA)=zi_transient_F(1,end,AA);
    plot_absorber_h(iteration,1,AA)=ab_transient_F+plots(iteration,1,AA);                           % time-dependent absorber position
        else
%     'start falling'
    
        % 속도가 클 때 뒤집기 위해서



        %% Forces
        %% Drag force

        for i=1:node_MATRA-1                                       %for all height of FAST ( mesh point)
            alpha=floor(zi_transient_F(1,i,AA)/dz)+1;                   %coolant node 의 정수
            gamma=zi_transient_F(1,i,AA)/dz;
            beta=floor(zi_transient_F(1,i,AA)/dz);                      %coolant node 의 정수
            if beta <=0
                beta=beta+1;
                alpha=alpha+1;
                gamma=gamma+1;
            end
            Shear_Drag(i,1)=(mu_c(beta)*(alpha-gamma)+mu_c(alpha)*(gamma-beta) ) *(-v_i_p(1,1,AA)+... %surface drag force 
                (v_i_p(2,1,AA)))/dr;
        end
        %% 나중에 온도차 있을 때 shear_Drag 다시 mesh마다
        %% 밀도가 다르므로 그것을고려해야 함. sum(shear)*dA 와 shear *A_side 는 많은 차이점을 보임.        
        %%shear stress
        F_drag(AA)    =sum(Shear_Drag(1:end,1))*dA;         
        %% pressure force
        P_drop_com(AA)=0.42*(((D)/(D_p))^2)  * den_c(1,1) * 1/2 * v(AA)^2 ;    %입구영역에서 compaction에 의한 압력강하
        dP(AA)    =Pe_dz(AA)*H - P_drop_com(AA)-(max(v_i_p(:,:,AA))-min(v_i_p(:,:,AA)))^2/2*den_c(1,1);                      %(윗면압력-밑면압력)에 작용하는 압력차이(dynamic pressure는 서로 같다고 가정)
        F_pressure(AA)= -dP(AA)*A_front;

                %% bouyancy
        for i=1:node_MATRA-1
                alpha=floor(zi_transient_F(1,i,AA)/dz)+1;                   %coolant node 의 정수
                gamma=zi_transient_F(1,i,AA)/dz;
                beta=floor(zi_transient_F(1,i,AA)/dz);                      %coolant node 의 정수
                if beta <=0
                    beta=beta+1;
                    alpha=alpha+1;
                    gamma=gamma+1;
                end
                if i==1
                    Fb(i,AA)=-(den_c(beta)*(alpha-gamma)+den_c(alpha)*(gamma-beta))*g*A_front*dh*0.5;
                end
                if i==node_MATRA-1
                    Fb(i,AA)=-(den_c(beta)*(alpha-gamma)+den_c(alpha)*(gamma-beta))*g*A_front*dh*0.5;
                end
                Fb(i,AA)=-(den_c(beta)*(alpha-gamma)+den_c(alpha)*(gamma-beta))*g*A_front*dh;
        end
        F_buoyancy(AA)=sum(Fb(:,AA));
        %gravity
        F_gravity =density_FAST* V * g ;
        %% velocity & acceleration 


         acceleration=(F_buoyancy(AA)+F_drag(AA)+F_gravity+F_pressure(AA))/density_FAST/ V; %가속도
         a(AA)=acceleration;                                                
         v(AA)=v(AA)+a(AA)*dt;                                                     %속도
         s(AA)=s(AA)+v(AA)*dt;                                                     %위치(적분된 위치)
         zi_transient_F(:,:,AA)=zi_F+(H_p-H)+s(AA);
        %iteration 횟수
        %% each time step data saving in transient
    plota(iteration,1,AA)=a(AA);
    plotv(iteration,1,AA)=v(AA);
    plots(iteration,1,AA)=s(AA);
    plotb(iteration,1,AA)=F_buoyancy(AA);
    plotg(iteration,1,AA)=F_gravity;
    plott(iteration,1,AA)=trigger(AA);
    plotd(iteration,1,AA)=F_drag(AA);
    plotp(iteration,1,AA)=F_pressure(AA);
    plotq(iteration,1,AA)=abs(v(AA))*A_front;                             %mass flow FAST pusing
    plotpd(iteration,1,AA)=Pe_dz(AA);
      %  zi_transient=zi_transient+s;
      %actual FAST position and each meshing point
    plot_reactivity(iteration,1)=total_reactivity;
    plot_feedback(iteration,1)=feedback_reactivity;
    plot_inserted_r(iteration,1)=inserted_reactivity_FAST;
    plot_position_start(iteration,1,AA)=zi_transient_F(1,1,AA);
    plot_position_end(iteration,1,AA)=zi_transient_F(1,end,AA);
    plot_absorber_h(iteration,1,AA)=ab_transient_F+plots(iteration,1,AA);                         % time-dependent absorber position
                       % time-dependent absorber position
end
        %역 방향으로 위치가 잡혔을 때 ( 위쪽으로 부딪힌다는 의미 ) 
    

if  zi_transient_F(1,end,AA) >= H_p                                          %벽에 다시 닿고 힘을 위로 받을 때
    %% each time step data saving in steady state
    a(AA)=0;
    v(AA)=0;
    s(AA)=0;
    F_drag(AA)=0;
    F_pressure(AA)=0;
    Q_error(AA)=0;
    Pe_dz(AA)=0;
    zi_transient_F(:,:,AA)=zi_F+(H_p-H)+s(AA);

        
    plota(iteration,1,AA)=a(AA);
    plotv(iteration,1,AA)=v(AA);
    plots(iteration,1,AA)=s(AA);
    plotb(iteration,1,AA)=F_buoyancy(AA);
    plotg(iteration,1,AA)=F_gravity;
    plott(iteration,1,AA)=trigger(AA);
    plotd(iteration,1,AA)=F_drag(AA);
    plotp(iteration,1,AA)=F_pressure(AA);
    plotq(iteration,1,AA)=abs(v(AA))*A_front;                             %mass flow FAST pusing
    plotpd(iteration,1,AA)=Pe_dz(AA);
      %  zi_transient=zi_transient+s;
      %actual FAST position and each meshing point
    plot_reactivity(iteration,1)=total_reactivity;
    plot_feedback(iteration,1)=feedback_reactivity;
    plot_inserted_r(iteration,1)=inserted_reactivity_FAST;
    plot_position_start(iteration,1,AA)=zi_transient_F(1,1,AA);
    plot_position_end(iteration,1,AA)=zi_transient_F(1,end,AA);
    plot_absorber_h(iteration,1,AA)=ab_transient_F+plots(iteration,1,AA);                         % time-dependent absorber position
end

if v(AA)>0
    reverse=-1;
else
    reverse=1;
end

    
    
%% Mass flow and Boundary condition change
%양수일 때도 고려하기 위해서


Q=-v(AA)*A_front;                                                          %volumetric flow(아래 밀어낸 값을 양으로 본다)

% 경계조건
v_i(1,1,AA)=v(AA);
v_i(end,1,AA)=0;
% 초기 속도장 계산
for i = 2 :num_gap-1
    v_i(i,1,AA)=-dt/den_c(1,1)*Pe_dz(AA) + mu_c(end,1)*dt/den_c(1,1)/dr^2*(ri(i)+dr/2)/ri(i)*v_i_p(i+1,1,AA) ...
        +  mu_c(end,1)*dt/den_c(1,1)/dr^2*(ri(i)-dr/2)/ri(i)*v_i_p(i-1,1,AA)...
        + (-2*mu_c(end,1)*dt/den_c(1,1)/dr^2+1)*v_i_p(i,1,AA);
end




for i=1:num_gap-1
    c1(i)=2*pi*ri(i)*dr*(v_i(i,1,AA)+...
        v_i(i+1,1,AA))/2;                         %sum of two velocity for differential flow area . page 39 reference
end
Q_e=sum(c1);                                                               %estimated flow(오른쪽으로 흐르고 있는 양)
Q_error=Q-Q_e;                                                             %오른쪽으로 흘러 가야 하는 양
Pe_dz_a=10;
Pe_dz_b=-10;
%% Pressure gradient correction technique
if Q == 0
    Q_error=0;
end
while abs((Q_error)/Q)>=tol2                                           %밀어내는 양과 옆으로 들어가는 양이 같을 때까지 P의 변화를 준다.
    %첫번 째 pressure gradient 에 대해서
    for i = 2 :num_gap-1
    v_i(i,1,AA)=-dt/den_c(1,1)*Pe_dz_a + mu_c(end,1)*dt/den_c(1,1)/dr^2*(ri(i)+dr/2)/ri(i)*v_i_p(i+1,1,AA) ...
        +  mu_c(end,1)*dt/den_c(1,1)/dr^2*(ri(i)-dr/2)/ri(i)*v_i_p(i-1,1,AA)...
        + (-2*mu_c(end,1)*dt/den_c(1,1)/dr^2+1)*v_i_p(i,1,AA);
    end

    for i=1:num_gap-1
        c1(i)=2*pi*ri(i)*dr*(v_i(i,1,AA)+...
            v_i(i+1,1,AA))/2;                         %sum of two velocity for differential flow area . page 39 reference
    end
    Q_e_a=sum(c1);                                                               %estimated flow(오른쪽으로 흐르고 있는 양)
    Q_error_a=Q-Q_e_a;                                                             %오른쪽으로 흘러 가야 하는 양
    RE_a=Q_error_a/Q;                
    %맞으면 탈출
    if abs(RE_a)<=tol2 
        Pe_dz(AA)=Pe_dz_a;

        break
    end
    
    %두번 째 pressure gradient 에 대해서
    
    for i = 2 :num_gap-1
    v_i(i,1,AA)=-dt/den_c(1,1)*Pe_dz_b + mu_c(end,1)*dt/den_c(1,1)/dr^2*(ri(i)+dr/2)/ri(i)*v_i_p(i+1,1,AA) ...
        +  mu_c(end,1)*dt/den_c(1,1)/dr^2*(ri(i)-dr/2)/ri(i)*v_i_p(i-1,1,AA)...
        + (-2*mu_c(end,1)*dt/den_c(1,1)/dr^2+1)*v_i_p(i,1,AA);
    end

    for i=1:num_gap-1
        c1(i)=2*pi*ri(i)*dr*(v_i(i,1,AA)+...
            v_i(i+1,1,AA))/2;                         %sum of two velocity for differential flow area . page 39 reference
    end
    Q_e_b=sum(c1);                                                               %estimated flow(오른쪽으로 흐르고 있는 양)
    Q_error_b=Q-Q_e_b;                                                             %오른쪽으로 흘러 가야 하는 양
    RE_b=Q_error_b/Q;         
    if abs(RE_b)<=tol2 
        Pe_dz(AA)=Pe_dz_b;
        break
    end
    
    % 둘 사이에 있는 경우, 
    if Q_error_b*Q_error_a<0                                               % 차이가 음수면 둘 사이에 있다고 판단.
        if abs(Q_error_a)>=abs(Q_error_b)
            Pe_dz_a=(Pe_dz_b+Pe_dz_a)/2;
        else
            Pe_dz_b=(Pe_dz_b+Pe_dz_a)/2;
        end
    end
        
        
    if Q_error_b*Q_error_a>0 
        if abs(Q_error_a)>=abs(Q_error_b)
            Pe_dz_a=(Pe_dz_b-Pe_dz_a)/(Q_e_b-Q_e_a)*(Q-Q_e_b)+Pe_dz_b;
        else 
            Pe_dz_b=(Pe_dz_b-Pe_dz_a)/(Q_e_b-Q_e_a)*(Q-Q_e_a)+Pe_dz_a;
        
        end
    end
    
    
    
        
    
           
end % pressure gradient correction method

%% 



v_i_p(:,:,AA)=v_i(:,:,AA);                                                                 %Previous velocity filed를 저장해 둠

% figure(2)
% showtime2=num2str(Ldt);                                   %round-off error removal. Fourth decimal place
% format short                                                 
% plot(ri,v_i(1:num_gap,1))
% title('Velocity profile')
% ylabel('Velocity[m/s]')
% ylim([-0.01 0.01]);
% xlabel('radius[m]')
% legend(strcat('time=',showtime2))
% drawnow;

end

if (-2*mu_c(end,1)*dt/den_c(1,1)/dr^2+1)<0
    'Reduce Time step '
    break
end

    
% store steady-state fuel and coolant temperature
Tf_mean_transient=mean(mean(mean(Plot_TC(1:node_f,:,:))));                        %fuel mean temperature at steady-state
Tc_mean_transient=mean(mean(mean(Plot_TC(node_t+1,:,:))));                        %coolatn mean temperature at steady-state
Ts_max_transient=max(max(Plot_TC(node_t+1,:,:)));                                %Coolant maximum temperature

% AA 그룹에서 의 온도 저장 
plot_Tfuel(iteration,1,AA)= Plot_TC(1,H_a/ddfh+1,AA);                          % Fuel temperature at Z=H_a
plot_Tfuel(iteration,2,AA)= Plot_TC(1,ceil(H_a/2/ddfh),AA)  ;                 
plot_Tfuel(iteration,3,AA)= Plot_TC(1,1,AA) ;                                 
plot_Tcladding(iteration,1,AA)=  Plot_TC(node_f,H_a/ddfh+1,AA)   ;                          
plot_Tcladding(iteration,2,AA)=  Plot_TC(node_f,ceil(H_a/2/ddfh),AA) ;                      
plot_Tcladding(iteration,3,AA)=  Plot_TC(node_f,1,AA)   ;                                   
plot_Tcoolant(iteration,1,AA)=Plot_TC(node_t+1,H_a/ddfh+1,AA);
plot_Tcoolant(iteration,2,AA)=Plot_TC(node_t+1,ceil(H_a/2/ddfh),AA);
plot_Tcoolant(iteration,3,AA)=Plot_TC(node_t+1,1,AA);
plot_Tpin(iteration,1,AA)= Plot_TC(node_t+1+node_r-floor((D_p/2-t_cf)/ddr),H_a/ddfh+1,AA)   ;                                             % pin  temperature Z = Ha (active core height), Z=Ha/2, Z= 0
plot_Tpin(iteration,2,AA)= Plot_TC(node_t+1+node_r-floor((D_p/2-t_cf)/ddr),ceil(H_a/2/ddfh),AA);                                                % pin  temperature Z = Ha (active core height), Z=Ha/2, Z= 0
plot_Tpin(iteration,3,AA)= Plot_TC(node_t+1+node_r-floor((D_p/2-t_cf)/ddr),1,AA);                                                % pin  temperature Z = Ha (active core height), Z=Ha/2, Z= 0
plot_Tpincoolant(iteration,1,AA)= Plot_TC(end,H_a/ddfh+1,AA);
plot_Tpincoolant(iteration,2,AA)= Plot_TC(end,ceil(H_a/2/ddfh),AA) ;           %floor((D_p/2-t_cf)/ddr)  is  the allocated node number to pin 
plot_Tpincoolant(iteration,3,AA)= Plot_TC(end,1,AA);      


%Power control by reactivity
%inserted reactivity by the FAST

%%  Inserted length calculation for inner FUal assembly
for AA  = 1 : length(Group_Q)
Absorber_Bottom = plot_absorber_h(iteration,1,AA);

absorber_section=Absorber_Bottom;
for i = 1 : absorber_section_number-1
    absorber_section= [absorber_section,Absorber_Bottom+h_ab/absorber_section_number*i]; 
end
% Inserted length calculation
Inserted_length=[];
for i = 1 : length(absorber_section)
    if H_actual-absorber_section(i)>0
        Temporary_Length=H_actual-absorber_section(i);
    end
    if H_actual-absorber_section(i)<=0
        Temporary_Length=0;
    end
    if H_actual-absorber_section(i)>=h_ab/absorber_section_number
        Temporary_Length=h_ab/absorber_section_number;
    end
        
    
    Inserted_length=[Inserted_length,Temporary_Length];
end
        
inserted_reactivity_FAST(AA)=Group_IC(AA)*sum(Absorber_Reactivity_Ratio.*Inserted_length);
end
inserted_reactivity_FAST=sum(inserted_reactivity_FAST);

%Rod withdrawal event
Reactivity_by_CWE=rod_withdrawal_reactivity/rod_withdrawal_time*Ldt; % 시간당 제어봉이 삽입하는 반응도량
if Ldt>=rod_withdrawal_time
    Reactivity_by_CWE=rod_withdrawal_reactivity;
end


%negative feedback by fuel and cladding
feedback_reactivity=(neg_f)*(Tf_mean_transient-Tf_mean_steady)+ (neg_c)*(Tc_mean_transient-Tc_mean_steady)...
    +(neg_a+neg_r+neg_CEDL)*(Tc_mean_transient-Tc_mean_steady);
%inserted total reactivity
total_reactivity=(inserted_reactivity_FAST)+(feedback_reactivity)+(Reactivity_by_CWE);
% fprintf('inserted_reactivity_FAST: %g [cent] ,\n  Absorber height : %g \n', inserted_reactivity_FAST,plot_absorber_h(iteration,1) )

% Boundary condition
MatrixPKE=zeros(7,7);
MatrixPKE(1,1)=1-(total_reactivity-PKE.beta_eff)/GT*dt;
for i = 1 : 6
    MatrixPKE(1+i,1+i)= 1+PKE.ramda(i)*dt;
    MatrixPKE(1+i,1)= -PKE.beta(i)/GT*dt;
    MatrixPKE(1,1+i)= -PKE.ramda(i)*dt;
end

% Implicit scheme
MatrixA=[Neutron(iteration);Precursors(iteration,1);Precursors(iteration,2);Precursors(iteration,3);Precursors(iteration,4);Precursors(iteration,5);Precursors(iteration,6)];
Next_Step=MatrixPKE\MatrixA;
Neutron(iteration+1)=Next_Step(1);
Precursors(iteration+1,1)= Next_Step(2);
Precursors(iteration+1,2)= Next_Step(3);
Precursors(iteration+1,3)= Next_Step(4);
Precursors(iteration+1,4)= Next_Step(5);
Precursors(iteration+1,5)= Next_Step(6);
Precursors(iteration+1,6)= Next_Step(7);



%Power change
q_t=q*transient_condition*Neutron(iteration+1).*Group_Q(:,:);
% BBB=G*A_c_f*2*(Plot_TC(node_t+1,1)*cp_coolant(1,1)-Plot_TC(node_t+1,end)*cp_coolant(end,1))
% fprintf('power: %g [W] ,\n  heat from coolant : %g [W]', q_t,AAA )


plot_Tem_f(iteration,1)=Tf_mean_transient;
plot_Tem_c(iteration,1)=Tc_mean_transient;
plot_Tem_cm(iteration,1)=max(max(Plot_TC(node_f+1,1:end,:)));
plot_Tem_fm(iteration,1)=max(max(Plot_TC(1,1:end,:)));
plot_Tem_s(iteration,1)=Ts_max_transient;

%% data saving session
average_channel_order= 97;                                                       %nomalized power 가 1 인 지점
aco=average_channel_order;

% 몫과 나머지
share=fix(iteration_main/iteration_cut);
rest=rem(iteration_main,iteration_cut);
if rest==0       %이터레이션 돌고, 한번 더돌았을 때, 
    share
    datestr(now)

    data_filename=['To' num2str(Ldt) 'Sec_final.mat'];
    save(data_filename)
    iteration=1;
    %이전 값을 이어 붙이기
    for AA = 1 : Group_N
    plota(iteration,1,AA)=plota(iteration_cut,1,AA);
    plotv(iteration,1,AA)=plotv(iteration_cut,1,AA);
    plots(iteration,1,AA)=plots(iteration_cut,1,AA);
    plotb(iteration,1,AA)=plotb(iteration_cut,1,AA);
    plotg(iteration,1,AA)=plotg(iteration_cut,1,AA);
    plott(iteration,1,AA)=plott(iteration_cut,1,AA);
    plotd(iteration,1,AA)=plotd(iteration_cut,1,AA);
    plotp(iteration,1,AA)=plotp(iteration_cut,1,AA);
    plotq(iteration,1,AA)=plotq(iteration_cut,1,AA);
    plotpd(iteration,1,AA)=plotpd(iteration_cut,1,AA);
    plot_reactivity(iteration,1)=plot_reactivity(iteration_cut,1);
    plot_feedback(iteration,1)=plot_feedback(iteration_cut,1);
    plot_inserted_r(iteration,1)=plot_inserted_r(iteration_cut,1);
    plot_position_start(iteration,1,AA)=plot_position_start(iteration_cut,1,AA);
    plot_position_end(iteration,1,AA)=plot_position_end(iteration_cut,1,AA);
    plot_absorber_h(iteration,1,AA)=plot_absorber_h(iteration_cut,1,AA);
    Neutron(iteration,1)=Neutron(iteration_cut,1);
    Precursors(iteration,:)=Precursors(iteration_cut,:);
    
    plot_den_b(iteration,AA)=plot_den_b(iteration_cut,1,AA);
    plot_den_t(iteration,AA)=plot_den_t(iteration_cut,1,AA);
    plot_xe(iteration,AA)=plot_xe(iteration_cut,1,AA);
    plot_Tfuel(iteration,1:3,AA)=plot_Tfuel(iteration_cut,1:3,AA);
    plot_Tcladding(iteration,1:3,AA)=   plot_Tcladding(iteration_cut,1:3,AA);
    plot_Tcoolant(iteration,1:3,AA)=plot_Tcoolant(iteration_cut,1:3,AA);
    plot_Tpin(iteration,1:3,AA)=plot_Tpin(iteration_cut,1:3,AA);
    plot_Tpincoolant(iteration,1:3,AA)=plot_Tpincoolant(iteration_cut,1:3,AA);
    
    end
    iteration=0;
    Ldt
end

iteration=iteration+1;
iteration_main=iteration_main+1;
Ldt=iteration_main*dt;

if floor(Ldt/0.01)==Ldt/0.01
    Ldt
    datestr(now)
    to=toc;
    (50-Ldt)*to/0.01/3600
    '시간남음'
    tic
end

end
