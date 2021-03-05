%% EXPLANATIONS 
% version v5(:v4) : Full undampD & nonlin(cubicForce) symbolic equations are included in A and B. 
% zeta=0 and gamma = 0 anyway for natural freq maps.
%
% The below is signed frequency calculation and plot for the linear, undamped , unforced overhung 
% rotor system. (see eg. Zilli's work)
% 
% 1) As Dr Shaw did in 2019 Normal Forms article, 
% firstly below,
% -the rotating frame eigenvalue prob is solved. Mode shapes and natFreqs are found.
% -The resulting responses are checked for whirl direction in rotating frame. 
% -Then the natFreqs are signed accordingly. 
% -Lastly, the stationary frame whirl speeds are calculated from 
%    wSta_calc = wRot_sign + Omeg. 
% -wRot_sign and |wSta_calc| (OPTIONally wSta_calc) are plotted.  
% 
% 2) Secondly, the same is done in stationary frame: 
% -the stationary frame eigenvalue prob is solved. Mode shapes and natFreqs are found.
% -The resulting responses are checked for whirl direction in stationary frame. 
% -Then the natFreqs are signed accordingly: wSta_signed
% -Lastly, the rotating coordinate whirls speeds are calculated form 
%    wRot_calc = wSta_sign - Omeg. 
% -|wSta_sign| (OPTIONally wSta_sign) and wRot_calc are plotted.
%
% THE SAME RESULTS ARE OBTAINED, going from either stationary or rotating frame systems.

%% ==================
%===ROTATING=FRAME===
%====================
clear,
syms Omeg Ip zeta gamma r2

n = 2 ; 
I = eye(2,2) ; Zero = zeros(2,2); J = [0,-1;1,0] ; 

A = [ -Omeg*J*(Ip-2)+2*zeta*I,   I
                I            ,  Zero ] ;
B = [ Omeg^2*(Ip-1)*I+I+2*zeta*Omeg*J+gamma*r2*I , Zero
                   Zero                          ,  -I  ] ;

zeta_val = 0.0 ; 
gamma_val = 0.0 ; 
A = subs(A,[zeta gamma],[zeta_val gamma_val]) ;
B = subs(B,[zeta gamma],[zeta_val gamma_val]) ;

[V,D] = eig(-A^(-1)*B) ; 
sRot = simplify(diag(D)) 

OmegRange = 0.1:0.1:7 ; 
for i = 1:length(OmegRange)
  Omeg_val = OmegRange(i) ;
  Ip_val = 0.143 ;

  sRot_val = double(subs(  sRot, [Omeg Ip], [Omeg_val Ip_val]  )) ;
  V_val    = double(subs(  V   , [Omeg Ip], [Omeg_val Ip_val]  )) ;
  %| Adjust the sRot order up to n, then adjust the modeShape matrix V accordingly.
  for j = 2:n
    if sRot_val(j)==conj(sRot_val(j-1)) 
      sRot_val([j-1+n,j]) = sRot_val([j,j-1+n]) ;
      V_val(:,[j-1+n,j]) = V_val(:,[j,j-1+n]) ;
      V_val([j-1+n,j],:) = V_val([j,j-1+n],:) ;
    end 
  end
  %| Adjust the sRot order from n+1 to 2*n, then adjust the modeShape matrix V accordingly.
  for j = n+1:2*n 
    if sRot_val(j) ~= sRot_val(j-n) 
      index = find(sRot_val==conj(sRot_val(j-n)));
      sRot_val([index,j]) = sRot_val([j,index]) ;
      V_val(:,[index,j]) = V_val(:,[j,index]) ;
      V_val([index,j],:) = V_val([j,index],:) ;
    end
  end                                               
%   wRot_val = abs(real(sRot_val))/zeta_val %|Damping present: s=-w*zeta+j*w_d , w_d = w*\/1-Zeta^2  
%   wRot_val = abs(imag(sRot_val)/sqrt(1-zeta_val^2)) 
  wRot_val = abs(sRot_val) ;%|No damping in natFreq map anyway.
  
  for r = 1:n
    T{r} = 2*pi./wRot_val ;
    t{r} = 0:0.01:T{r} ;
    
    part{1} = V_val(1:2,r  )*exp(sRot_val(r  )*t{r}) ;
    part{2} = V_val(1:2,r+n)*exp(sRot_val(r+n)*t{r}) ;
    respRot{r} = part{1} + part{2} ;
    
    angleStart = atan2( respRot{r}(2,1),respRot{r}(1,1) ) + 2*pi*(respRot{r}(2,1)<0) ;
    angleNext  = atan2( respRot{r}(2,5),respRot{r}(1,5) ) + 2*pi*(respRot{r}(2,5)<0) ;
    angDiff = angleNext-angleStart ;
    if angDiff>pi || angDiff<0 
      sign(r) = -1 ;
    else 
      sign(r) = +1 ;
    end
  end
  signMat = diag(sign) ;
  
  wRot_val_sign(:,i) = signMat*wRot_val(1:2) ;
  wSta_calc(:,i) = wRot_val_sign(:,i)+OmegRange(i) ;
  wSta_calc_killSign(:,i) = abs(wSta_calc(:,i)) ;
end

figure
ax1 = subplot(121) ; hold on ; grid on 
ax2 = subplot(122) ; hold on ; grid on
%| OPTION
%| :1
% plot(ax1,OmegRange, wSta_calc(1,:), "b", ...
%          OmegRange, wSta_calc(2,:), "r")
%| :2
plot(ax1,OmegRange, wSta_calc_killSign(1,:), "b", ...
         OmegRange, wSta_calc_killSign(2,:), "r")
%| :.
plot(ax2,OmegRange, 1*wRot_val_sign(1,:), "b-" , ...
         OmegRange, 2*wRot_val_sign(1,:), "b--", ...
         OmegRange, 3*wRot_val_sign(1,:), "b:" , ...
         OmegRange, 4*wRot_val_sign(1,:), "b-.", ...
         OmegRange, 1*wRot_val_sign(2,:), "r-" , ...
         OmegRange, 2*wRot_val_sign(2,:), "r--", ...
         OmegRange, 3*wRot_val_sign(2,:), "r:" , ...
         OmegRange, 4*wRot_val_sign(2,:), "r-.")

plot(ax1,OmegRange,OmegRange)
title(ax2,["Rot"," "])
axis(ax1,"square")
axis(ax2,"square")

xlabel(ax1,"$$\it \hat{\Omega} $$","interpreter","latex","fontsize",14)
xlabel(ax2,"$$\it \hat{\Omega} $$","interpreter","latex","fontsize",14)

ylabel(ax1,"$$\it \hat{\omega_r} $$","interpreter","latex","fontsize",14)
ylabel(ax2,"$$\it \hat{\omega_r} $$","interpreter","latex","fontsize",14)

legend(ax1,["Forward","Backward","1st EO"])
legend(ax2,["Forward","2xForward","3xForward","4xForward", ...
            "Backward","2xBackward","3xBackward","4xBackward"])



%% =================================================================================================
%===STATIONARY=FRAME================================================================================
%===================================================================================================
%|Nt: More Commenting here.

clear,clc
syms Omeg Ip zeta gamma r2

n = 2 ;
I = eye(2,2) ; Zero = zeros(2,2); J = [0,-1;1,0] ; 

%| The X-EOM system's linear(gamma=0),unforced(eH=0),undamped(zeta=0) state-space form: Ax+B=0
A = [-Omeg*Ip*J+2*zeta*I,  I
              I         , Zero ] ;
B = [ I+gamma*r2*I, Zero
          Zero    ,  -I ] ;
zeta_val = 0.0 ; 
gamma_val = 0.0 ; 
A = subs(A,[zeta gamma],[zeta_val gamma_val]) ;
B = subs(B,[zeta gamma],[zeta_val gamma_val]) ;
   
%| EIGENVAL PROBLEM : s*A*x + B*x = 0  |=>  s*A*V = -B*V*D 
%| ...where D = [s1 0;0 s2] and V = [{x1} {x2}] 
   
for i = 1
%|___AAA> METHOD-A : EIGENVAL PROB : [s*A+B]*x=0 
% C = (A*s + B) ;
% detC = det(C) ;
% eq = detC == 0 ; 
% s_vec = solve(eq,s) ;
% %| Now, eigenvalues appear as below. (I just did manually 1i=sqrt(-1) ) ;
% s_vec = [-1i*(1+(Ip^2*Omeg^2)/2+(Ip*Omeg*(Ip^2*Omeg^2+ 4)^(1/2))/2)^(1/2)
%          -1i*(1+(Ip^2*Omeg^2)/2-(Ip*Omeg*(Ip^2*Omeg^2+ 4)^(1/2))/2)^(1/2)
%           1i*(1+(Ip^2*Omeg^2)/2+(Ip*Omeg*(Ip^2*Omeg^2+ 4)^(1/2))/2)^(1/2)
%           1i*(1+(Ip^2*Omeg^2)/2-(Ip*Omeg*(Ip^2*Omeg^2+ 4)^(1/2))/2)^(1/2) ];
% w_vec = s_vec*1i ;
% w_vec12 = w_vec(1:2) ;%|Unsigned frequencies 
% % V ... ... ...
%|___AAA.
end  %|Method A

%|___BBB> METHOD-B : EIGENVAL PROB : s*x = [-A^(-1)*B]*x  |=>  [-A^(-1)*B]*V = V*D
[V,D] = eig(-A^(-1)*B) ;
sSta = diag(D) ;
%|___BBB.

%|___YYY> 
OmegRange = 0.1:0.1:7 ; 
for i = 1:length(OmegRange)
  Omeg_val = OmegRange(i) ;
  Ip_val = 0.143 ;
  
  sSta_val = double(vpa(subs(  sSta, [Omeg Ip], [Omeg_val Ip_val]  ))) ;%|Subs>calcs in sym>calcs.
  V_val    = double(vpa(subs(  V   , [Omeg Ip], [Omeg_val Ip_val]  ))) ;%|Subs>calcs in sym>calcs.
  %| Adjust the sSta order up to n, then adjust the modeShape matrix V accordingly: [[x],[x]*]
  for j = 2:n %|MkSr dt if drs a conj pair in 1st n el of s, replace one by 
    if sSta_val(j)==conj(sSta_val(j-1)) 
      sSta_val([j-1+n,j]) = sSta_val([j,j-1+n]) ;
      V_val(:,[j-1+n,j]) = V_val(:,[j,j-1+n]) ;
      V_val([j-1+n,j],:) = V_val([j,j-1+n],:) ;
    end 
  end
  %| Adjust the sSta order from n+1 to 2*n, then adjust the modeShape matrix V accordingly.
  for j = n+1:2*n 
    if sSta_val(j) ~= sSta_val(j-n) 
      index = find(sSta_val==conj(sSta_val(j-n))); %|No tolerance needed.
      sSta_val([index,j]) = sSta_val([j,index]) ;
      V_val(:,[index,j]) = V_val(:,[j,index]) ;
      V_val([index,j],:) = V_val([j,index],:) ;
    end
  end
%   wSta_val = abs(real(sSta_val))/zeta_val %|Damping present: s=-w*zeta+j*w_d , w_d = w*\/1-Zeta^2 
%   wSta_val = abs(imag(sSta_val))/sqrt(1-zeta_val^2) ;
  wSta_val = abs(sSta_val) ;%|No damping in natFreq map anyway.

  %| Find sign of the Omeg_val in hand
  for r = 1:n %|mode
    T{r} = 2*pi./abs(wSta_val(r)) ;
    t{r} = 0:0.01:T{r} ;
    
    %| Calc response in the Omeg_val in hand
    part{1} = V_val(1:2,r  )*exp(sSta_val(r  )*t{r}) ;
    part{2} = V_val(1:2,r+n)*exp(sSta_val(r+n)*t{r}) ;
    respSta{r} = part{1} + part{2} ; 

    %| Obserb in which direction d traj is rotatJ | add 2pi if atan2 returns a negative angle.
    angleStart = atan2( respSta{r}(2,1),respSta{r}(1,1) ) + 2*pi*(respSta{r}(2,1)<0) ;
    angleNext  = atan2( respSta{r}(2,5),respSta{r}(1,5) ) + 2*pi*(respSta{r}(2,5)<0) ;
    angDiff = angleNext - angleStart ;
    if  angDiff>pi || angDiff<0
      sign(r) = -1 ;
    else 
      sign(r) = +1 ;
    end
  end
  signMat = diag(sign) ;
  
  wSta_val_sign(:,i) = signMat * wSta_val(1:2) ;%|Signed natFreqs | Alrdy 2dof, so manually 1:(n=2)
  wRot_calc(:,i) = wSta_val_sign(:,i) - OmegRange(i) ;
  wSta_val_killSign(:,i) = abs(wSta_val_sign(:,i)) ;%| wSta_val_killSign for traditional Campbell
end 

%| All for figure below: 
figure
ax1 = subplot(121) ; hold on 
ax2 = subplot(122) ; hold on
%| OPTION 
%| :1: Non traditional Campbell (signed sta plot) 
plot(ax1,OmegRange, wSta_val_sign(1,:), "b", ...
         OmegRange, wSta_val_sign(2,:), "r")
%| :2: Traditional Campbell plot (no sign in sta plot)
% plot(ax1,OmegRange, wSta_val_killSign(1,:), "b", ...
%          OmegRange, wSta_val_killSign(2,:), "r")
%| :.

%|
plot(ax2,OmegRange, 1*wRot_calc(1,:), "b-" , ...
         OmegRange, 2*wRot_calc(1,:), "b--", ...
         OmegRange, 3*wRot_calc(1,:), "b:" , ...
         OmegRange, 4*wRot_calc(1,:), "b-.", ...
         OmegRange, 1*wRot_calc(2,:), "r-" , ...
         OmegRange, 2*wRot_calc(2,:), "r--", ...
         OmegRange, 3*wRot_calc(2,:), "r:" , ...
         OmegRange, 4*wRot_calc(2,:), "r-.")
plot(ax1,OmegRange,OmegRange)
title(ax1,["Sta"," "])
axis(ax1,"square")
axis(ax2,"square")
%|YYY.


%|BOTTOM.