clc
clear all
close all
% Import data--------------------------------------------------------------
adrs='E:\Research_for_ThesiS\ATF_Simulation_Code\IPIN_Competition\IPIN 2017\DataSim';
adr2='E:\Research_for_ThesiS\ATF_Simulation_Code\IPIN_Competition\IPIN 2017\Code_Sim\CAR\Validation\Madgwick_Comp';
load([adrs, '/Wall_m']);%save('Wall_m')
addpath([adr2,'/common_func/']);addpath([adr2,'/Madgwick/']);
RGB_paint = imread('CAR_paint_control.jpg');load('controlpoint');%
%----- pixel to meter -----------------------------------------------------
scale=0.05207600; lonX= -3.48367648;latY= 40.31320308;angle=-8.77680000;[Xcent_m,Ycent_m]=ll2utm(latY,lonX);
[X_cp Y_cp]=pixel2Meter(controlpoint,RGB_paint,angle,scale,Xcent_m,Ycent_m);
clear RGB_paint;RGB_paint = imread('CAR_paint_refined.jpg');%imshow(RGB_paint)
%controlpoint(65,:)=[860 288]; save('controlpoint')
%% ----------------   lat lon 2 meter  -------------------------------------
DataLog=load([adrs, '/logfile_CAR_R02-2017_A5']);
%% ----------------   lat lon 2 meter  -------------------------------------
[P0x P0y]=ll2utm([40.31320308 -3.48367648]);angle=-8.77680000;scale=0.05207600;% meter/pixel    % lon degree
[x,y]=ll2utm(DataLog.Posi(:,3),DataLog.Posi(:,4));
if size(DataLog.Posi,1)==22
    y(8,1)=4.462634891186919e+06;y(11,1)=4.462634966491276e+06;y(15,1)=4.462646467531692e+06;
end
%%----------------- for each path between two ref. points -----------------
start=15
finish=16;%estmated_path_LV(:,3)=Pdr_new(3,:)';%A5_3_13t14=estmated_path_LV';
%A5_2_15t16=Pdr_new;
%save('A5_2_15t16')
%--------------------- OLD PDR estmated     -------------------------------
old_pdr=PDR_system_A5(x,y,start,finish,DataLog.Acce,DataLog.Posi,DataLog.Gyro,DataLog.Magn,DataLog.Ahrs,DataLog.Wifi,DataLog.Ble4,DataLog.Ligh);
%---------------------  New PDR estmated    -------------------------------
[Pdr_new, OLD,SL_teta]=PDR_new(old_pdr,finish,x,y);PDR_NEW_orginal=Pdr_new;PDR_old_orginal=OLD;
%SL_teta(1,1)=SL_teta(1,2);SL_teta(2,1)=SL_teta(2,2);
% ---------------- Error of Landmark point --------------------------------
ex=x(finish)- old_pdr.coor(end,1);ey=y(finish)- old_pdr.coor(end,2);clear ex ey
%% ========================================================================
plot(Wall_m(:,1),Wall_m(:,2),'.'); hold on;plot(Pdr_new(1,:),Pdr_new(2,:),'k.-','MarkerSize',4);hold on
plot(OLD(1,:),OLD(2,:),'b.-','MarkerSize',4);hold on;plot(x,y,'o:','MarkerSize',3);hold on
for i=1:size(x,1)
    text(x(i,1),y(i,1),num2str(i),'FontSize',8,'Color','k');
end;hold on;plot(X_cp(:,1), Y_cp(:,1),'^')
%-------- Start Particle Filter -------------------------------------------
alpha=0.2;% fusion coef  >>>>  a(D)+(1-a)d
var=.5;
np=50;% particle number
M0=10; % number of selected importance particle
mu =[Pdr_new(1,1) Pdr_new(2,1)];      % Gaussian point generation
Particle=Gaussian_propag(mu,var,np);
Particle=Inithial_corrected_particle(Particle,mu,RGB_paint,Wall_m);% corrected invalid particle
%hold on;plot(Particle(:,1),Particle(:,2),'.')
clear mu 
%% ----------------------   PF Loop  --------------------------------------
for stp_id=1:size(old_pdr.Sl,2) % for all step detected on current path
    stp_id %stp_id=stp_id+1   
    Prtcl_new=prtcl_loc_update(Particle,SL_teta(1,stp_id),SL_teta(2,stp_id));
    %hold on;plot(Prtcl_new(:,1),Prtcl_new(:,2),'.')
    %hold on;plot(Particle(:,1),Particle(:,2),'.')
%-----------  on wall all test ---------------------------------------------
    On_wall = Wall_on_all(Prtcl_new,RGB_paint);
    onw_id0 = find(On_wall==1); 
%     if size()
%-----------  wall crossing test -------------------------------------------
    for jj=1:size(Prtcl_new,1) 
        flg(jj)=Wall_Crossing(Prtcl_new(jj,:),Particle(jj,:),RGB_paint);
    end;  id_cr=find(flg==1);
    sum_mask=On_wall+flg;
    invalid_particle(stp_id)=size(find(sum_mask~=0),2);clear sum_mask 
% --------------- Observation Error ----------------------------------------
er_invalid_State(stp_id)=invalid_particle(stp_id)/size(Prtcl_new,1);
if er_invalid_State(stp_id)>=0.8
    Error_count=[stp_id er_invalid_State(stp_id)]
    er_invalid_State(stp_id);
    s_Loc=Pdr_new(:,1);
    if stp_id==1
        past_loc=Pdr_new(1:2,1);
    else
        past_loc=Est_LVar(stp_id-1,1:2);
    end
    ns=1;%hold on;plot(Prtcl_new(:,1),Prtcl_new(:,2),'.')
    if er_invalid_State(stp_id)>0.9
        P=Gaussian_propag(past_loc,1,np);
        P=Inithial_corrected_particle(P,past_loc,RGB_paint,Wall_m);
        Prtcl_new=P;
    else
        P=Prtcl_new;
    end
    [cond_pos,New_loc_past]=Observation_Error(s_Loc,[X_cp Y_cp],past_loc',P,ns);clear P
%     if(Wall_Crossing(cond_pos,New_loc_past,RGB_paint)==0)
%         New_loc_past=cond_pos;
%     end
    if Wall_on_all(New_loc_past,RGB_paint)==1
        New_loc_past=cond_pos;
    end
    hold on;plot(cond_pos(1,1),cond_pos(1,2),'k*')
    hold on;plot(New_loc_past(1,1),New_loc_past(1,2),'go');clear past_loc
    if stp_id==1
        before_coef=0; 
    else
        before_coef=1;
    end
    Est_LVar(stp_id-1,1:2)=New_loc_past;
    [new_rout,Ltheta]=Corected_Path2(New_loc_past,stp_id,Pdr_new,before_coef);clear New_loc
    Est_LVar(1:stp_id-1,1:2)=new_rout(1:2,2:stp_id)';% past correction
    OLD=Pdr_new;
    Pdr_new=new_rout;  % Pdr_new correction
    hold on;plot(Pdr_new(1,:),Pdr_new(2,:),'r.-')
    SL_teta=Ltheta;
    Mu=Est_LVar(stp_id-1,1:2);
    Prtcl_new=Gaussian_propag(Mu,var,np);
    Prtcl_new=Inithial_corrected_particle(Prtcl_new,Mu,RGB_paint,Wall_m);
    %hold on;plot(Prtcl_new(:,1),Prtcl_new(:,2),'.')
    clear Mu New_rout Ltheta ns past_loc s_Loc    
else
    %--------------wall crossing corection -------------------------------------
    if size(id_cr,2)~=0
        wcr_loc=Prtcl_new(id_cr' ,:); 
        summask=On_wall+flg;vid=find(summask==0);
        if size(vid,2)~=0
            valid_loc=Prtcl_new(vid' ,:);
            pos_cr=Cross_2_closest(wcr_loc,valid_loc);
            Prtcl_new(id_cr' ,:)=pos_cr; 
            clear wcr_loc pos_cr vid summask flg id_cr
        end
    end
% ---------- wall on correction-------------------------------------------
    if size(onw_id0,2)~=0
        onwall_prt=Prtcl_new(onw_id0' ,:);nowall_id0=find(On_wall==0);nowall=Prtcl_new(nowall_id0' ,:);
        loc0=OnWall2closestpoint(onwall_prt,nowall);
        Prtcl_new(onw_id0' ,:)=loc0;
    end ;clear onw_id0 nowall_id0 loc0 onwall_prt  On_wall
   %hold on;plot(Prtcl_new(:,1),Prtcl_new(:,2),'g.') 
% ------------- distance correction ---------------------------------------
    idex_wdis=Dis_Part2Wallpoint(Prtcl_new,Wall_m,0.25);
    %hold on; plot(Prtcl_new(:,1),Prtcl_new(:,2),'g.')  
    if size(idex_wdis,1)~=0
        mask=zeros(size(Prtcl_new,1),1);mask(idex_wdis)=1;mask=mask.*(1:size(Prtcl_new,1))';
        walclos_loc=Prtcl_new(idex_wdis,:);valid_id=find(mask==0);
        valid_loc=Prtcl_new(valid_id,:);
        if size(valid_loc,1)~=0
            Pos_cls=Closewall_2_closest(walclos_loc,valid_loc) ;
            Prtcl_new(idex_wdis,:)=Pos_cls;
            clear walclos_loc idex_wdis  valid_loc mask Pos_cls
        end  %hold on; plot(Prtcl_new(:,1),Prtcl_new(:,2),'m.')
    end% end wall on correction
end %  if (Error_number(stp_id)/size(Prtcl_new,1))>=0.8
%--------------- Wheiting Step -------------------------------------------- 

    W=Weghting_prt(Prtcl_new,OLD(1:2,stp_id+1)',...
        Pdr_new(1:2,stp_id+1)',size(Prtcl_new,1),alpha,old_pdr.Sl(1,stp_id));
%----- 1th without resampling >>> sum W*Loc -------------------------------
    Est_W(stp_id,1:2)=sum(W.*Prtcl_new,1);
 %----- 2th importance methode --------------------------------------------
    [~,iDx]=sort(W,'descend');           % decressing
    M=min(M0,size(Prtcl_new,1));
    import_prtcl=Prtcl_new(iDx(1:M),:);
    Est_Import(stp_id,1:2)=[sum(import_prtcl(:,1))/M,sum(import_prtcl(:,2))/M];
    clear import_prtcl    
%----- 3th resampling and Estimated Loc.-----------------------------------
    index_rs = lowVarianceRS(1:size(Prtcl_new,1) , W, size(Prtcl_new,1));
    p_lv=Prtcl_new;%clear Est_W
    p_lv=p_lv(index_rs',:);
    %Prtcl_new=p_lv;
    Est_LVar(stp_id,1:2)=(sum(p_lv,1))/size(p_lv,1); 
    LV(stp_id,1:2)=Est_LVar(stp_id,1:2);
    %hold on; plot(Est_LVar(stp_id,1),Est_LVar(stp_id,2),'ko')
    clear index_rs p_lv

%----- Error Detection from PF estimated ----------------------------------
    
    Dif_W(stp_id,1)=sqrt(sum((Est_W(stp_id,1:2)-Pdr_new(1:2,stp_id+1)').^2,2));
    Dif_LV(stp_id,1)=sqrt(sum((Est_LVar(stp_id,1:2)-Pdr_new(1:2,stp_id+1)').^2,2));
    Dif_imp(stp_id,1)=sqrt(sum((Est_Import(stp_id,1:2)-Pdr_new(1:2,stp_id+1)').^2,2));
    er_tr=0.2;
    if Dif_LV(stp_id,1)>=er_tr % dif between pf est and pdr_new
%----------% wall on cheking of curent estimated --------------------------
        FLG_wo=Wall_on(Est_LVar(stp_id,1:2),RGB_paint);
        if (stp_id==1 && FLG_wo==1) % 1 >> invalid Location --- new particle for correction
            Mu0 =[Est_LVar(stp_id,1) Est_LVar(stp_id,2)];
            %Mu0 =[Pdr_new(1,1) Pdr_new(2,1)];
            Particle0=Gaussian_propag(Mu0,var,np);
            onwall_p=Wall_on_all(Particle0,RGB_paint);vwal=find(onwall_p==1);
            Particle0(vwal',:)=[];
            sum_dif2=sum((Particle0 - repmat(Mu0,size(Particle0,1),1)).^2,2);
            disss=sqrt(sum_dif2);[nm,idx_t]=min(disss);
            Est_LVar(stp_id,:)=Particle0(idx_t,:);% corrected is exchanged
            Mu =[Est_LVar(stp_id,1) Est_LVar(stp_id,2)];
        elseif FLG_wo==0
            Mu =[Est_LVar(stp_id,1) Est_LVar(stp_id,2)];
        end
         Mu =[Est_LVar(stp_id,1) Est_LVar(stp_id,2)];
         before_coef=1;
%         if stp_id==1
%             before_coef=1; % only curent step should be corrected 
%         else
%             before_coef=1; % all  pdr_new should be corrected 
%         end
        %hold on; plot(Est_LVar(stp_id,1),Est_LVar(stp_id,2),'ko')
        [new_rout,Ltheta]=Corected_Path2(Mu,stp_id,Pdr_new,before_coef);
        Est_LVar(1:stp_id,1:2)=new_rout(1:2,2:stp_id+1)';% past correction
        OLD=Pdr_new;
        Pdr_new=new_rout;  % Pdr_new correction
        hold on;plot(Pdr_new(1,:),Pdr_new(2,:),'r.-')
        SL_teta=Ltheta;
       Prtcl_new=Inithial_corrected_particle(Prtcl_new,Mu,RGB_paint,Wall_m);
       clear Mu Ltheta  new_rout FLG_wo
       %hold on;plot(Prtcl_new(:,1),Prtcl_new(:,2),'.')
    end % tr
%----------------- Error cheking of two point estimated -------------------
    if stp_id>1   
% Error % wall crossing cheking of two estimated point
% ??????????????????????????????????? new path az divar na pf estimation
%         ErrEst(stp_id,1:2)=Error_alarm(Est_LVar(stp_id-1,1:2),...
%             Est_LVar(stp_id,1:2),RGB_paint,stp_id); 
          ErrEst(stp_id,1:2)=Error_alarm(Pdr_new(1:2,stp_id)',...
              Pdr_new(1:2,stp_id+1)',RGB_paint,stp_id);
        if ErrEst(stp_id,1)~=0  % Error detected
            ns=1;
            %[cond_loc, est_cont_LV]=Corected_Loc([X_cp Y_cp],Est_LVar(stp_id,1:2),Prtcl_new,1);
            Mu=Est_LVar(stp_id-1,1:2);%hold on;plot(Est_LVar(stp_id,1),Est_LVar(stp_id,2),'k.')
            %Mu=Pdr_new(1:2,stp_id)';
% !!!!!!!!!!!!!!!
            sigma=SL_teta(2,stp_id);
            Prtcl_new=Gaussian_propag(Mu,sigma,np);
            %hold on;plot(cond_loc(:,1),cond_loc(:,2),'c*')
            %Prtcl_new=Inithial_corrected_particle(Prtcl_new,Mu,RGB_paint,Wall_m);clear Mu
            [cond_loc, est_cont_LV]=Corected_Loc2(Pdr_new(:,1),SL_teta,[X_cp Y_cp],Est_LVar(stp_id,1:2)',Prtcl_new,ns);
%----------- % LOS between cond_loc , est_cont_LV  --------------------------
            if Wall_Crossing(cond_loc,est_cont_LV,RGB_paint)==1
                Prt=Gaussian_propag(est_cont_LV,var,np);
                for jj=1:size(Prt,1)
                    flg(jj)=Wall_Crossing(Prt(jj,:),cond_loc,RGB_paint);
                end;  id_cr=find(flg==1);
                if size(id_cr,2)~=0
                    valid_L=Prt(id_cr,:);
                    dist2past=sqrt(sum( (valid_L- repmat(cond_loc,size(valid_L,1) ,1)).^2  ,2) );
                    [~,min_id]=min(dist2past);
                    new_P=valid_L(min_id,:);
                    est_cont_LV=new_P;clear new_P min_id dist2past valid_L flg id_cr
                end
            end
%----------- %End  LOS between cond_loc , est_cont_LV  --------------------------
%             hold on;plot(est_cont_LV(1,1),est_cont_LV(1,2),'g*')
            coef=1;
            [new_Path lt]=Corected_Path2(est_cont_LV,stp_id,Pdr_new,coef);
            Est_LVar(1:stp_id,1:2)=new_Path(1:2,2:stp_id+1)';
            Mu =[est_cont_LV(1,1) est_cont_LV(1,2)];
            OLD=Pdr_new;
            Pdr_new=new_Path;  % Pdr_new correction
            hold on;plot(cond_loc(1,1),cond_loc(1,2),'g*','MarkerSize',5)
            hold on;plot(Pdr_new(1,:),Pdr_new(2,:),'g.-');
            SL_teta=lt;
            Prtcl_new=Gaussian_propag(Mu,var,np);
            Prtcl_new=Inithial_corrected_particle(Prtcl_new,Mu,RGB_paint,Wall_m);
            % corrected invalid particle  
            clear Mu lt new_Path est_cont_LV
        end  
    end   
%----------------------------- PLOT -----------------------------------
    %plot(Est_Import(stp_id,1 ),Est_Import(stp_id, 2),'g*:') ;hold on
      plot(Est_LVar(stp_id,1),Est_LVar(stp_id,2),'ks:');hold on
%     %plot(Est_W(stp_id,1),Est_W(stp_id,2),'co:');hold on
      drawnow
%----------------- New Particle exchanging -------------------------------
    Particle=Prtcl_new;clear Prtcl_new
end 
%% ====================== PLOT ===========================================
figure(2)
estmated_path_LV(1,1:2)=Pdr_new(1:2,1);% start point
estmated_path_LV(2:size(LV,1)+1,1:2)=LV(:,1:2);
plot(Wall_m(:,1),Wall_m(:,2),'.'); hold on
% estmated_pathLV(1,1:2)=Pdr_new(1:2,1);% start point
% estmated_pathLV(2:size(Est_LVar,1)+1,1:2)=Est_LVar;
%hold on;plot(estmated_pathLV(:,1),estmated_pathLV(:,2),'g.:')
hold on;plot(Pdr_new(1,:),Pdr_new(2,:),'k.:');hold on;
%plot(PDR_NEW_orginal(1,:),PDR_NEW_orginal(2,:),'r.:');hold on;
plot(PDR_old_orginal(1,:),PDR_old_orginal(2,:),'r.:')
%legend('map','reconstructed Path','raw Path')
hold on;plot(estmated_path_LV(:,1),estmated_path_LV(:,2),'g.:');hold on
plot(x,y,'o:','MarkerSize',3);hold on
for i=1:size(x,1)
    text(x(i,1),y(i,1),num2str(i),'FontSize',8,'Color','k');
end
%legend('map','reconstructed Path','raw Path','ref.')


%% Error  ------------------------------------------------
% id_refr=start+1;%clear cur_loc_ref ;clear dt;clear  diff_2
% for io=1:finish-start-1
%     %ii=ii+1;
%     ref_time=DataLog.Posi(id_refr,1)
%     cur_loc_ref(io,1:2)=[x(id_refr,1) y(id_refr,1)];
%     est_time=Pdr_new(3,2:end)';
%     dt=abs(repmat(ref_time,size( est_time,1) , 1)-est_time);
%     [val id_srt]=sort(dt);
%     ref_estima(io,1:2)=Pdr_new(1:2,id_srt(1)+1);
%     diff_2=ref_estima(io,:)- cur_loc_ref(io,1:2);
%     error(io,1)=sqrt(sum((diff_2).^2))
%     id_refr=id_refr+1;
%     clear val id_srt dt ref_time  est_time 
% end
% plot(ref_estima(:,1),ref_estima(:,2),'b+')
% hold on 
% plot(cur_loc_ref(:,1),cur_loc_ref(:,2),'go')
% legend('map','reconstructed Path','raw Path','ref.','estimated','test.point')
% clear  ref_estima















