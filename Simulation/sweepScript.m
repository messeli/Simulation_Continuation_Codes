clc, clear,% close all
clockONplot = tic ;

%% Solve for a range of Omeg (rotor speed)

%|___DEF___|> Define solver and coords
solvers = ["mySolver", "ode45"] ;
frames = ["sta" "rot"] ;
%|:stationary {phixH; phixHP; phiyH; phiyHP}
%|:rotating   {    u;      v; u_dot;   vdot}
solver = solvers(2) ;
%|___DEF___|. 

%|Nt:The initialInitial condition qn_ to start the sineSweep in phixH phiyH phixHP phiyHP coords.
%|...The regular initialInitialCond is [0.6;0;0;0.5] 
% qn_=[-1.1379;-1.9816;-7.7839;5.3304];%|singleLoop sweep contin, obD from qn_=[0;0.6;-1.246;0] @4.81
qn_ =[0.3808;-1.1155;1.6958;2.6338] ;%|doubleLoop sweep contin, obD from qn_=[0;0.6;-1.246;0] @2.91

tend = 3000 ;%4000
tol = 1e-7 ;%1e-9
Omeg_range_txt = "2.91:0.01:2.93" ;%"4.81:-0.01:0.00";%"0.0:0.1:7" ;%"5.60:0.01:5.65"
Omeg_range = eval(Omeg_range_txt) ;
isBack = false ;%|ACTION:Backward sine sweep toggle. False means forward.

iRange = 1:length(Omeg_range) ;
if isBack 
  iRange = flip(iRange) ;
end

for frame = frames(2)
  frame
  qn = qn_ ;%|Set initCond of the first Omeg run as the initialInitialCond
  
  for i = iRange
    %|___INIT___|> INITIAL CONDITION
    Omeg = Omeg_range(i) ;
    if frame == "sta"
      %| ACTION:If phi_EOM is UNCOMMENTED in func_...
      "pass";
      %| ACTION:If X_EOM is UNCOMMENTED in func_...
      % T_phi2X = [0 -1 0 0;1 0 0 0;0 0 0 -1;0 0 1 0];%|q_phi = T_phi2X . q_X
      % qn = T_phi2X^(-1) * qn ;%|q_X = T_phi2X^(-1) . q_phi
    elseif frame == "rot"
      T_X2rot = [1 0 0 0;0 1 0 0;0 -Omeg 1 0;Omeg 0 0 1] ;%|q_X = T_X2rot . q_rot
      T_phi2X = [0 -1 0 0;1 0 0 0;0 0 0 -1;0 0 1 0] ;%|q_phi = T_phi2X . q_X
      qn = T_X2rot^(-1) * T_phi2X^(-1) * qn; %|q_u = T_X2rot^(-1) . T_phi2X^(-1) . q_phi
    end
    %|___INIT___|.    

    %|___SOL___|> SELECT SOLVER - SOLVE
    switch solver
      case "mySolver" %|Still not quite functional 
      [T{i},Q{i},XisContact{i},hMat{i}] = func_mysolver(Omeg,qn,tend,tol,frame) ;
      
      case "ode45"
      [T{i},Q{i}] = func_ode45(Omeg,qn,tend,tol,frame) ;
      hMat{i} = [];
      for j = 1:length(T{i})-1
        hMat{i}(j) = T{i}(j+1)-T{i}(j) ; 
      end
    end
    
    if frame == "sta"
      % qn = T_phi2X * Q{i}(:,end) ;%|ACTION:Function in X-EOM | X to phi 
      %|:SWEEP:return the last phi-coord, it will take care for X-frame init cond
      qn = Q{i}(:,end) ;%|ACTION:Function in phi, so return in phi-coords.
    elseif frame == "rot"
      qn = T_phi2X * T_X2rot * Q{i}(:,end) ;%|rot to phi
      %|:SWEEP:return the last phi-coord, it will take care for rot-frame init cond
      % qn = qn + 0.5*max(abs(qn))*(2*rand(4,1)-1) ;%|Random noise(!) at the next Omeg init cond by 1%
      % qn = qn*1.29;
    end
    %|___SOL___|.

  end %|Entire Omeg_range finished for one coordinate system

  %| SIGNAL AT FINISH OF CURRENT COORD
  clockOFFplot = toc(clockONplot) 
  disp(['COMPLETED in ', num2str( clockOFFplot / 60 ), ' mins.'])


  
  %|___FIG___|> PLOTTING CURRENT COORD RESULTS
  %|___FIG_INDIV___|> INDIVIDUAL PLOT 
  for Omeg_check = [2.91,2.92] %Omeg_range(1)
    I = find( abs( Omeg_range-Omeg_check ) <= 0.001 );
    individualplot_ISOcubicStiffnessNoContact(T{I},Q{I},...
                                            Omeg_check,[0.80 1.0],false); 
    %|:TF - Animation
  end
  %|___FIG_INDIV___|.

  %|___FIG_3D2D___|> 3D&2D PLOTS, POINCARE SECTIONS
  clockONplot = tic ; 
  fprintf("\nPlotting...") 
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
    from = round( 0.95 * length(Q{i}) , 0 ) ;% 0.97
    to = length(Q{i}) ;
    Omeg = Omeg_range(i) ;

    %|___FIG_3D2D_3D___|> 3D PLOT - STACK
    plot3( ax1, Q{i}(1,from:end), Q{i}(2,from:end), Omeg*ones(1,length(from:to)), ...
           "k" ) %"k."
    %|___FIG_3D2D_3D___|.

    %|___FIG_3D2D_2D___|> 2D PLOT - AMPLITUDE
    rH = sqrt( Q{i}(1,:).^2 + Q{i}(2,:).^2 ) ;
    plot( ax2, Omeg*ones(1,length(from:to)) , rH(from:to) , "k.", ...
               Omeg, max(rH(from:to)), "bo")
    %|___FIG_3D2D_2D___|.

    %|___FIG_3D2D_POIN___|> STROBOSCOPIC POINCARE SECTION 
    %|___FIG_3D2D_POIN_INDX___|> GRAB PERIOD-CLOSE INDEXES 
    TPeriod = 2*pi/Omeg ;
    k = 1 ; r = 1 ; IPeriod = [] ; 
    fromPoin = round( 0.60 * length(Q{i}) , 0 ) ;% 0.97
    toPoin = length(Q{i}) ;    
    while k*TPeriod <= T{i}(end)
      if k*TPeriod > T{i}(fromPoin) && k*TPeriod <= T{i}(toPoin)
        [~,IPeriod(r)] = min( abs( T{i} - k*TPeriod ) ) ;
        r = r+1 ;%|IPeriod index
      end
      k = k+1 ;%|Period's number.
    end
    %|___FIG_3D2D_POIN_INDX___|.
    
    plot3( ax1, Q{i}(1,IPeriod), Q{i}(2,IPeriod), Omeg*ones(1,length(IPeriod)),...
           ".","lineWidth",1,"color",colour) % ".","lineWidth",2,"color","k")
    plot ( ax2 , Omeg*ones(1,length(IPeriod)) , rH(IPeriod) , ...
           ".","lineWidth",1,"color",colour) %"bo"
    %|___FIG_3D2D_POIN___|. POINCARE
  end %|Continue to the next Omeg in Omeg_range
  %|___FIG_3D2D___|. 3D&2D 
  clockOFFplot = toc(clockONplot) 
  fprintf("Plotted in %.1f min", clockOFFplot/60)
  %|___FIG___|. PLOTTING
  
end  %|Switch to the next coord frame
