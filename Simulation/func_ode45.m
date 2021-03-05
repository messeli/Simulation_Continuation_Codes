function [T,Q] = func_ode45(Omeg,qn,tend,tol,coordFrame) 
% Omeg = 5 ;
% qn = [0.6;0;0;0.5] ;
% tend = 3000 ;
% tol = 1e-7 ;% That is RelTol of ode45.
% coordFrame= 'sta' ;% otherwise 'rot'

if nargin <= 4  %|number of arguments input
  f = @f_sta ;%|Default coordinate.
else
  if     coordFrame == "sta", f = @f_sta ;
  elseif coordFrame == "rot", f = @f_rot ;
  end
end

phi = 0      ;
hTry  = 1e-8   ;
tn    = 0      ; 
tnm1  = 0      ;
t     = [tn ]  ;
q     = [qn' ] ;
T     = [t ]   ; 
Q     = [q ]   ;
gamma = 0.25   ;%0.25 %|f3=gamma*f, f3=k3*cC^3, f=kr*cC, cC-clearance of ZilliSystem. 
%|:gamma is direct meas f k3, if cC assumd cnst.
mH    = 0.9    ;%0.9
isInitialStep = true ; 

tic
disp(['STARTED  Omeg = ', num2str(Omeg)])

%|___AAA___|> CALL ODE45
options = odeset('RelTol',tol,'AbsTol',tol*1e-2,'initialStep',hTry, ...
                                                'events',@shootAwayEvent) ;
[t,q] = ode45(f,[tn tend],qn,options) ;

%| Try to eliminate the following transpose for future versions. DEBUG
Q = q  ;  T = t   ;% Each row is a data point, so long cols.
Q = Q.';  T = T.' ;% Each col is a data point, so long rows. 
%|___AAA___|. 

TOC = toc ;
disp(['FINISHED Omeg = ', num2str(Omeg), ' in ',num2str(TOC),' sec.'])


%% %%%%%%%%%-%%
%%-FUNCTIONS-%%
%%-%%%%%%%%%-%%
%% STA
function dq = f_sta(tn,qn)
  % fprintf("%6.2f of %4.2f\n",tn/tend*100,Omeg) % Track progress.
  %| cC is clearance, but c can be cos. So, it is cC not to confuse.
  
    %| Prepare Omeg and tnm1
  if ~isInitialStep  % If not the initial step. 
    hnm1 = (tn-tnm1)  ;% h for prev step used to calc Omeg for this step.
    phi = phi + hnm1*Omeg  ;% tnm1: time atW prev step was calcd.
    tnm1 = tn  ;%|next step preparation. 
  else 
    phi = phi + (hTry/4)*Omeg ; 
    %|:The specified 1st step size is always devided by 4, why is that??, see in hMat??
    isInitialStep = false ;
  end
  
  zetas = [0.010,0.008,0.005,1e-3,1e-4,1e-5] ;
  JpH = 0.143 ; epsH = 3.53e-1 ; zeta = zetas(1) ; OmegP = 0 ; 

  %| Trasnform btw d two stationary frames: X_hat = ang_y_hat ; Y_hat = -ang_x_hat ;
  staTypes = ["ang_EOM", "X_EOM"] ;
  %| ang_EOM {ang_x_hat; ang_y_hat; ang_x_hat_dot; ang_y_hat_dot}.
  %| X_EOM {X_hat, Y_hat, X_hat_dot, Y_hat_dot}.
  staType = staTypes(2) ;%|ACTION 
  switch staType
    case "ang_EOM"
    q1 = qn(1) ; q2 = qn(2) ; q3 = qn(3) ; q4 = qn(4) ;%|Leave ang as ang
    case "X_EOM"
    q1 =-qn(2) ; q2 = qn(1) ; q3 =-qn(4) ; q4 = qn(3) ;%|Transform from X to ang
  end
  %|:ACTION:DONT FORGET to UNCOMMENT the relevant X_EOM initial conditions from Zilli_script_...m
  %|...if you chose to use X_EOM, vice versa.

  %|___XXX___|> The below are the equations of ang_EOM. 
  r2 = q1^2 + q2^2 ;

  B = +1 ;%|Backward Integration sign: -1 means backward integration is attempted. FAIL yet.
  dq1 = B*(B*q3 ) ;
  dq2 = B*(B*q4 ) ;
  dq3 = (- JpH * (B*Omeg) * B*q4 ...
        - 2 * (B*zeta) * (B*q3) ...
        - q1 ... 
        + mH * epsH * ( OmegP * cos(phi) - (B*Omeg)^2 * sin(phi) ) ...
        - gamma * r2 * q1 ) ;%|ISO Stiffness one.
  dq4 = (+ JpH * (B*Omeg) * (B*q3) ...
        - 2 * (B*zeta) * (B*q4) ...
        - q2 ...
        + mH * epsH * ( OmegP * sin(phi) + (B*Omeg)^2 * cos(phi) ) ...
        - gamma * r2 * q2 ) ;%|ISO Stiffness one.
  %|___XXX___|. ang_EOM

  switch staType 
    case "ang_EOM"
    dq = [dq1; dq2;dq3; dq4] ;%|Leave ang as ang
    case "X_EOM"
    dq = [dq2;-dq1;dq4;-dq3] ;%|Transform from ang to X
  end
end
%% ROT
function dq = f_rot(tn,qn)
  %|cC is clearance, but c can be cos. So, it is cC not to confuse.
  % fprintf("%6.2f of %4.2f\n",tn/tend*100,Omeg) % Track progress.

  %|Prepare Omeg and tnm1
  if ~isInitialStep  %|If not the initial step. 
    hnm1 = (tn-tnm1)  ;%|h for prev step used to calc time at this step.
    phi = phi + hnm1*Omeg  ;%|tnm1: time atW prev step was calcd.
    tnm1 = tn  ;%|next step preparation. 
  else %|If it is the initial step, update the angle with specified hTry. 
    phi = phi + (hTry/4)*Omeg ; 
    %|:The specified 1st step size is always devided by 4, why is that??, See in hMat
    isInitialStep = false ;
  end
  
  zetaDepot = [0.010,0.008,0.005,1e-3,1e-4,1e-5] ;
  JpH = 0.143 ; epsH = 3.53e-1 ; zeta = zetaDepot(1) ; OmegP = 0 ; 
  
  q1 = qn(1)     ; q2 = qn(2) ; q3 = qn(3)   ; q4 = qn(4) ;
  % u:rot x-axis ; u_prime    ; v:rot y-axis ; v_prime    ;
  
  %|___UUU___|> The below are the equations of ang_EOM. 
  r2 = q1^2 + q2^2 ;
  
  B = +1 ;%|backwardIntegrationSign: -1 means backward integration is attempted. FAIL yet.
  dq1 = B*(B*q3) ; 
  dq2 = B*(B*q4) ;    
  dq3 = (- (B*Omeg) * (JpH-2) * (B*q4) ...
        - 2 * (B*zeta) * (B*q3) ...
        + ( (B*Omeg)^2 * (1-JpH)  -  1 ) * q1 ...
        + 2 * (B*zeta) * (B*Omeg) * q2 ... %|This is - @2019 Shaw.
        + mH*epsH*(B*Omeg)^2 ...
        - gamma * r2 * q1 ) ;%|ISO Stiffness one.
  %     - gamma/4*( 3*q1*r2 + c4*q1*(r2-4*q2^2) + s4*q2*(r2-4*q1^2) ) ;
  %| Depot: 
  %   s = sin(phi) ; c = cos(phi) ;
  %   s4 = sin(4*phi) ; c4 = cos(4*phi) ;
  %| Simplified original independent stiffness one:
  %         - gamma/4*( (3+c4)*q1^3 + (s4)*q2^3 ) ...
  %         - gamma/4*( (-3*s4)*q1^2*q2 + (3-3*c4)*q1*q2^2) ;
  %| Original independent stiffness one:
  %         - gamma*( (c^4+s^4)*q1^3 + (c^3*s-c*s^3)*q2^3 ) ...
  %         - gamma*( (-3*c^3*s+3*c*s^3)*q1^2*q2 + (6*c^2*s^2)*q1*q2^2) ;
  dq4 = (+ (B*Omeg) * (JpH-2) * (B*q3) ...
        - 2 * (B*zeta) * (B*q4) ...
        + ( (B*Omeg)^2 * (1-JpH)  -  1 ) * q2 ...
        - 2 * (B*zeta) * (B*Omeg) * q1 ... %|This is + @2019 Shaw. 
        - mH*epsH*OmegP ...
        - gamma * r2 * q2 ) ;%|ISO Stiffness one.
  %     - gamma/4*( 3*q2*r2 + c4*q2*(r2-4*q1^2) - s4*q1*(r2-4*q2^2) ) ;

  %| Depot: 
  %| Simplified original independent stiffness one: 
  %         - gamma/4*( (3+c4)*q2^3 - (c^3*s-c*s^3)*q1^3 ) ...
  %         - gamma*(-(-3*s4)*q1*q2^2 + (3-3*c4)*q1^2*q2) ;
  %| Original independent stiffness one: 
  %         - gamma*( (c^4+s^4)*q2^3 - (c^3*s-c*s^3)*q1^3 ) ...
  %         - gamma*(-(-3*c^3*s+3*c*s^3)*q1*q2^2 + (6*c^2*s^2)*q1^2*q2) ;
  %|___UUU___|.
  dq = [dq1;dq2;dq3;dq4];
end

%% contactChange event
% 
% function [value,isterminal,direction] = myEventFcn(tn,qn)
%     value = sqrt( qn(1)^2 + qn(3)^2 ) - 1 ;
% %   isterminal = 0 ;% Not terminated as isContact is updated inside f.
%     isterminal = 1 ;% Terminated upon a contactChange to update isContact. 
%     direction = 0 ;
% end
%% ShootAway event
function [value,isterminal,direction] = shootAwayEvent(tn,qn)
  value = max(qn(1),qn(3)) - 1e3 ;%|Coordinates exceeds 100.0 normalised position. 
  isterminal = 1 ; 
  direction = 0;
end
    
%% myoutputfunc  !!! NOT UP-TO-DATE !!! 
% 
% function status = myoutputfunc(t,q,flag)
%     persistent hMat
%     switch flag
%         case 'init'
%             hMat = [hMat h];% Initially the h must be hDefault, otherwise it is the last h before contactChange
%         case ''
%             hMat = [hMat, t(end)-t(end(-1) ];
%         case 'done' % when it's done
%             assignin('base','hMat',hMat); % get the data to the workspace.
%     end
%     status = 0;
% end

end


