clc;clear all;close all
%% Add Pathes and Load Data ------------------------------------------------
adrs='.\DataSim';
adr2='...\functions';
addpath([adr2,'/common_func/']);addpath([adr2,'/Madgwick/']);addpath([adr2,'/MatData/']);
addpath([adr2,'/MapImages/']);
load([adrs, '/Wall_m']);load('PathN_pxm');load('controlpoint');
load('MatrixGraph');load('PathSs');load('controlpoint');
RGB_paint = imread('CAR_paint_control.jpg');
Data=load([adrs, '/logfile_CAR_R01-2017_S4']);
%% Ref. Coordinate lat., lon. to be transformed to meter  ---------------------------------------
ref_lat=Data.Posi(:,3);ref_lon=Data.Posi(:,4);[x_ref,y_ref]=ll2utm(ref_lat,ref_lon);
clear ref_lat ref_lon
%% pixel to meter of control points (which identify door and cross way locations on map -------------
scale=0.05207600; lonX_cnt= -3.48367648;latY_cnt= 40.31320308;angle=-8.77680000;
[Xcent_m,Ycent_m]=ll2utm(latY_cnt,lonX_cnt);
[X_cp Y_cp]=pixel2Meter(controlpoint,RGB_paint,angle,scale,Xcent_m,Ycent_m);
clear RGB_paint;RGB_paint = imread('CAR_paint_refined.jpg');
%% PDR for each path between two ref. points --------------------------
grayy=rgb2gray(RGB_paint);
start=1   % start refrence point 
finish=5  % end refrence point 
old_pdr=PDR_system(x_ref,y_ref,start,finish,Data.Acce,...
    Data.Posi,Data.Gyro,Data.Magn,Data.Ahrs,Data.Wifi,Data.Ble4,Data.Ligh);
%% New PDR estmated -------------------------------------------------------
[Pdr_new, OLD,SL_teta]=PDR_new(old_pdr,finish,x_ref,y_ref);
PDR_NEW_orginal=Pdr_new;
PDR_old_orginal=OLD;
%% Plot Wall,Ref., Old Pdr and New Pdr
plot(Wall_m(:,1),Wall_m(:,2),'.','MarkerSize',4);hold on;
plot(x_ref,y_ref,'o:','MarkerSize',3);hold on
plot(PDR_old_orginal(1,:),PDR_old_orginal(2,:),'m.-','MarkerSize',4);hold on
plot(Pdr_new(1,:),Pdr_new(2,:),'k.-','MarkerSize',4);hold on
for i=1:size(x_ref,1)
    text(x_ref(i,1),y_ref(i,1),num2str(i),'FontSize',8,'Color','k');
end;hold on;plot(X_cp(:,1), Y_cp(:,1),'g*')

%% Start Particle Filter --------------------------------------------------
% (0) initializing
alpha=0.4;% fusion coef  >>>>  a(D)+(1-a)d
var=.5;
np=50;    % particle number
M0=10;    % number of selected importance particle
mu =[Pdr_new(1,1) Pdr_new(2,1)];      % Gaussian point generation
Particle=Gaussian_propag(mu,var,np);
Particle=Inithial_corrected_particle(Particle,mu,RGB_paint,Wall_m);% corrected invalid particle
%hold on;plot(Particle(:,1),Particle(:,2),'.')
clear mu
%% ----------------------   PF Loop  --------------------------------------
for stp_id=1:size(old_pdr.Sl,2) % for all step detected on current path
    stp_id %stp_id=stp_id+1
% -------------- new location based observation ---------------------------
    Prtcl_new=prtcl_loc_update(Particle,SL_teta(1,stp_id),SL_teta(2,stp_id));
%------------- On Wall Detection ------------------------------------------
    On_wall = Wall_on_all(Prtcl_new,RGB_paint);
    onw_id0 = find(On_wall==1);
%-----------  Wall crossing test ------------------------------------------
    for jj=1:size(Prtcl_new,1)
        flg(jj)=Wall_Crossing(Prtcl_new(jj,:),Particle(jj,:),RGB_paint);
    end;  id_cr=find(flg==1);
    sum_mask=On_wall+flg;
    invalid_particle(stp_id,1)=size(find(sum_mask~=0),2);
    clear sum_mask
% ----------  Observation Error ( wall crossing + on wall) ----------------
    er_invalid_State(stp_id,1)=invalid_particle(stp_id,1)/size(Prtcl_new,1);
% ----------  Observation Error Correction --------------------------------
% if er_invalid_State(stp_id,1)>0.8
%     
% end
% ---------- wall on and wall crossing correction -------------------------      
if size(id_cr,2)~=0
    wcr_loc=Prtcl_new(id_cr' ,:);
    summask=On_wall+flg;
    vid=find(summask==0);
    if size(vid,2)~=0
        alrm_wc(stp_id,1)=1;
        valid_loc=Prtcl_new(vid' ,:);
        pos_cr=Cross_2_closest(wcr_loc,valid_loc);
        Prtcl_new(id_cr' ,:)=pos_cr;
        clear wcr_loc pos_cr vid summask flg id_cr
    end
end
% ---------- wall on correction--------------------------------------------
if size(onw_id0,2)~=0
    alrm_won(stp_id,1)=1;
    onwall_prt=Prtcl_new(onw_id0' ,:);nowall_id0=find(On_wall==0);
    nowall=Prtcl_new(nowall_id0' ,:);
    loc0=OnWall2closestpoint(onwall_prt,nowall);
    Prtcl_new(onw_id0' ,:)=loc0;
end ;clear onw_id0 nowall_id0 loc0 onwall_prt  On_wall
%--------------- Weighting Step -------------------------------------------
    W=Weghting_prt(Prtcl_new,OLD(1:2,stp_id+1)',...
        Pdr_new(1:2,stp_id+1)',size(Prtcl_new,1),alpha,old_pdr.Sl(1,stp_id));
%----- Low variance resampling and Estimated Loc --------------------------
    index_rs = lowVarianceRS(1:size(Prtcl_new,1) , W, size(Prtcl_new,1));
    p_lv=Prtcl_new;%clear Est_W
    p_lv=p_lv(index_rs',:);
    Prtcl_new=p_lv;
    Est_LVar(stp_id,1:2)=(sum(p_lv,1))/size(p_lv,1);
    LV(stp_id,1:2)=Est_LVar(stp_id,1:2);
    clear index_rs p_lv
    Particle=Prtcl_new;clear Prtcl_new
    drawnow
    hold on
    plot(LV(stp_id,1),LV(stp_id,2),'m.')
end
