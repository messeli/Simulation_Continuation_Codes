clear  % clc, close all

%|___DEF___|> DEFINE SOLVER AND COORDS
solvers = ["mySolver", "ode45"] ;
frames = ["sta" "rot"] ;
%|:stationary {phixH; phiyH; phixHP; phiyHP}
%|:stationary {   xH;    yH;    xHP;    yHP}
%|:rotating   {    u;      v; u_dot;   vdot}
solver = solvers(2) ;%|ACTION:Chose solver 1:mySolver 2:ode45
%|___DEF___|.


%|___INPUT___|> DEFINE SOLVER INPUTS
%     for randInitCondGenerationTimes = 1:5 %|ACTION:Do you need random init cond generation?
%       qn_ = rand(4,1)*20-10 ;%|+-20,15,10 done. 
% qn_ = [0.6;0;0;0.5] ;%|The default init cond I used most often is [0.6;0;0;0.5] 
qn_ = [0;0.6;-1.2460;0];%|The phi_initCond dt yields DOUBLE_LOOP at 2.91, SINGLE_LOOP at 4.81.
% qn_ = [-4.7626  ;   -3.2929   ;   3.5946  ;   -7.2689] ;
% qn_ = [0.79705      0.4986     -3.5505      3.5303]' ;
%|Nt:Define qn inside parfors sdt dy c handle it freely.
% qn_ = [0.3105;-7.2711;-1.7714;-1.9818] ;%|4.3's end point from tend=3000
tend = 3000 ;%|Switch tdir setting - FAIL:damping, See Mathworks question.
tol = 1e-7; 
Omeg_range_txt = "6.01:0.5:7.01" ;%"4.81" "0.01:0.1:7.01" ;%|Cld use "0.01:0.01:7.01" |"2.91"
Omeg_range = eval(Omeg_range_txt) ;
%|Nt: .dat files were generated at 4.05 for AUTO starting, for various zeta.
%|Nt: An UPO(unstable periodic orbit) was visualised at 7.01, See below Ctrl+F:UPO.
%|___INPUT___|.


for frame = frames(2) %|ACTION:Chose frames (1):sta (2):rot, NONE:both
  clockON = tic ;

  %|___LOOP___|> LOOP ALL Omeg_range FOR CURRENT COORD
  frame
  T = {};Q = {};
  p = 1 ; nIdxLeft = length(Omeg_range) ;

  while  nIdxLeft > 0
    %| Generate index batch for parfor
    iRange = p:p+min(8,nIdxLeft)-1 ; 
    disp(['NEW INDEX BATCH ' num2str(iRange)])

    parfor i = iRange  % parfor

      %|___LOOP_INIT___|> INITIAL CONDITION
      Omeg = Omeg_range(i) ;
      qn = qn_ ;%|Defined in phixH phixHP phiyH phiyHP coords.
      %|:You need to define qn here sdt each parfor can change it.
      if frame == "sta"
        %| ACTION:If phi_EOM is UNCOMMENTED in Zilli_func_...
        % "pass";
        %| ACTION:If X_EOM is UNCOMMENTED in Zilli_func_...
        T_phi2X = [0 -1 0 0;1 0 0 0;0 0 0 -1;0 0 1 0];%|q_phi = T_phi2X . q_X
        qn = T_phi2X^(-1) * qn ;
      elseif frame == "rot"
        T_X2rot = [1 0 0 0;0 1 0 0;0 -Omeg 1 0;Omeg 0 0 1] ;%|q_X = T_X2rot . q_rot
        T_phi2X = [0 -1 0 0;1 0 0 0;0 0 0 -1;0 0 1 0] ;%|q_phi = T_phi2X . q_X
        qn = T_X2rot^(-1) * T_phi2X^(-1) * qn; %|q_u = T_X2rot^(-1) . T_phi2X^(-1) . q_phi
      end
      %|___LOOP_INIT___|.


      %|___LOOP_SOL___|> SELECT SOLVER - SOLVE
      switch solver
        case "mySolver"
        [T{i},Q{i},XisContact{i},hMat{i}] = func_mysolver(Omeg,qn,tend,tol,frame) ;  

        case "ode45"
        Omeg = Omeg_range(i) ; 
        
        %|___LOOP_UPO> Visualise UPO>SPO |{u,v,du,dv}
        %| Case 1:A point on UPO at 7.01 in AUTO<<4.05's .dat contin w/ zeta=0.01 
        % qn = [-1;3.3567;22.2465;-1.46324] ;
        %| Case 2:A point on UPO at 3.47 in AUTO<<4.05's .dat contin w/ zeta=0.01
        % qn=[-0.64820093074;-8.6742629190E-2;-0.10752121904;0.79258623141];%|UPO correspJ t label29
        % qn=[-0.71363715809;-0.16475278048;-0.24866314284;1.07435451054];%|SPO correspJ t label31
        %|___LOOP_UPO.

        %qn=[-0.380009000620;-0.001297860010;0.000022150721;0.000206860717];%|{u,v,du,dv}
        %|:Dss lowAmp sol datum at Omeg=7.0, used fr AUTO IPS=1 "givenBackward" continuation. 
        [T{i},Q{i}] = func_ode45(Omeg,qn,tend,tol,frame) ;
        hMat{i} = [];
        for j = 1:length(T{i})-1
          hMat{i}(j) = T{i}(j+1)-T{i}(j) ; 
        end 
      end 
      %|___LOOP_SOL___|. 

    end %|This batch of indexRange finished.
    %| Prepare next indexRange for parfor |Check if parfor now works without this indexRange
    nIdxLeft = length(Omeg_range) - iRange(end); 
    p = iRange(end) + 1 ; 

  end %|Entire Omeg_range finished.
  %|___LOOP___|. FINISHED

  %| SIGNAL AT FINISH OF CURRENT COORD.
  clockOFF = toc(clockON) 
  disp([frame, 'COMPLETED in ', num2str( clockOFF / 60 ), ' mins.'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %|___FIG___|> PLOTTING CURRENT COORD RESULTS
  %|___FIG_INDIV___|> INDIVIDUAL PLOT 
  for Omeg_check = 6.51 %Omeg_range(1)
    I = find( abs( Omeg_range-Omeg_check ) <= 0.001 );
    individualplot_ISOcubicStiffnessNoContact(T{I},Q{I},...
                                            Omeg_check,[0.80 1.0],false); 
    %|:TF - Animation
  end
  %|___FIG_INDIV___|.

  %|___FIG_3D2D___|> 3D&2D PLOTS, POINCARE SECTIONS
  %|___FIG_3D2D_PREP___|> PREPARE FIGURE
  colour = [1.0 0.7 0.0];
  figure
  ax1 = subplot(121) ; hold(ax1,'on') ; grid(ax1,'on') ; axis(ax1, 'square','equal')
  ax2 = subplot(122) ; hold(ax2,'on') ; grid(ax2,'on') ; axis(ax2, 'square')
  
  if frame == "sta"
    xlabel(ax1,"$\it\hat{x}\,,\hat{\phi_y}$","interpreter","latex","FontSize",14)
    ylabel(ax1,"$\it\hat{y}\,,\hat{-\phi_x}$","interpreter","latex","FontSize",14) 
  elseif frame == "rot"
    xlabel(ax1,"$\it\hat{u}$","interpreter","latex","FontSize",14) 
    ylabel(ax1,"$\it\hat{v}$","interpreter","latex","FontSize",14)
  end
  zlabel(ax1,"Rotor speed","fontSize",10) %|Omega^H^a^t is Omeg
  title (ax1,[solver+" "+frame,"Omeg-range "+Omeg_range_txt," "],"FontSize",10)
  xlim(ax1,[-3 3]); ylim(ax1,[-3 3])
  %xlabel(ax2,"\Omega",'Interpreter','tex') %|Omega^H^a^t is Omeg
  xlabel(ax2,"$$\rm Rotor\,Speed,\,\it\hat{\Omega}$$",'Interpreter','latex',"FontSize",14) %|Omega^H^a^t=Omeg
  ylabel(ax2,"$\rm Amplitude,\,\it\hat{r}$","interpreter","latex","FontSize",14) 
  title(ax2,["phi-EOM initial conditions",num2str(qn_')," "],"fontSize",10)
  ylim(ax2,[0 5])
  %|___FIG_3D2D_PREP___|. 

  for i=1:length(Omeg_range)
    %| FIND INDEXES TO PLOT FROM PERCENTAGES
    from = round( 0.9959 * length(Q{i}) , 0 ) ;% 0.97
    to = length(Q{i}) ;
    Omeg = Omeg_range(i) ;

    %|___FIG_3D2D_3D___|> 3D PLOT - STACK
    plot3( ax1, Q{i}(1,from:end), Q{i}(2,from:end), Omeg*ones(1,length(from:to)), ...
           "k" ) %"k."
    %|___FIG_3D2D_3D___|.

    %|___FIG_3D2D_2D___|> 2D PLOT - AMPLITUDE
    rH = sqrt( Q{i}(1,:).^2 + Q{i}(2,:).^2 ) ;
    plot( ax2, Omeg*ones(1,length(from:to)) , rH(from:to) , "k.", ...
               Omeg , max(rH(from:to)), "bo") 
    %|___FIG_3D2D_2D___|.

    %|___FIG_3D2D_POIN___|> STROBOSCOPIC POINCARE SECTION 
    %|___FIG_3D2D_POIN_INDX___|> GRAB PERIOD-CLOSE INDEXES 
    TPeriod = 2*pi/Omeg ;
    k = 1 ; r = 1 ; IPeriod = [] ; 
    fromP = round( 0.60 * length(Q{i}) , 0 ) ;% 0.97
    toP = length(Q{i}) ;    
    while k*TPeriod <= T{i}(end)
      if k*TPeriod > T{i}(fromP) && k*TPeriod <= T{i}(toP)
        [~,IPeriod(r)] = min( abs( T{i} - k*TPeriod ) ) ;
        r = r+1 ;%|IPeriod index
      end
      k = k+1 ;%|Period's number.
    end
    %|___FIG_3D2D_POIN_INDX___|.
    
    plot3( ax1, Q{i}(1,IPeriod), Q{i}(2,IPeriod), Omeg*ones(1,length(IPeriod)),...
           ".","lineWidth",1,"color",colour) % ".","lineWidth",2,"color","k")
    plot(  ax2 , Omeg*ones(1,length(IPeriod)) , rH(IPeriod) , ...
           ".","lineWidth",1,"color",colour) %"bo"
    %|___FIG_3D2D_POIN___|. POINCARE
  end %|Continue to the next Omeg in Omeg_range
  %|___FIG_3D2D___|. 3D&2D 
  %|___FIG___|. PLOTTING
  
end %|Cycle to next coord system


%    end %|Cycle to new random initial condition