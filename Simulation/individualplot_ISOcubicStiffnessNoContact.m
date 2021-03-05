function individualplot_ISOcubicStiffnessNoContact(t,q,thetaP_check,fromTo_percentages, animateTF)
  % [from_percent to_percent] = fromTo_percentages 

  from = round( fromTo_percentages(1)*length(q), 0 ) ;if from==0, from=1; end
  to   = round( fromTo_percentages(2)*length(q), 0 ) ;

  %|___PREP> Prepare the figure
  figure , hold on 
  % viscircles([0 0],1,'color','b') ;
  plot(0,0,"r+")
  axis square
  xlim([-2.5 2.5]), ylim([-2.5 2.5]), xlabel("q_1"), ylabel("q_2")
  title(num2str(thetaP_check))
  plot(q(1,1),q(2,1),"ro","lineWidth",0.5) % Init point marked on plot. 
  %|___PREP. 

  %|___TRAJ> Plot the trajectory
  plot(q(1,from:to),q(2,from:to),"k") % X_Hat = phi_y ,  Y_Hat = -phi_x
  axis([-5 5 -5 5])
  %|___TRAJ.

  %|___POIN> Stroboscopic Poincare section 
  TPeriod = 2*pi/thetaP_check;
  k = 1 ; r = 1 ; IPeriod = [];
  while k*TPeriod <= t(end)
      if k*TPeriod > t(from) && k*TPeriod <= t(to)
          [~,IPeriod(r)] = min( abs(t-k*TPeriod) ) ;
          r = r+1 ;
      end
      k = k+1 ;
  end
  plot(q(1,IPeriod),q(2,IPeriod),"gs","lineWidth",2)
  %|___POIN. 

  %|___ANIM> Draw points to see the trajectory moving. 
  if animateTF == true 
      for i=from:to
          pause(0.01)
          if ( sqrt( q(1,i)^2 + q(2,i)^2 ) - 1  <  0 )
                  plot(q(1,i),q(2,i),"k.")
          else  
                  plot(q(1,i),q(2,i),"r.")
          end
      end
  end
  %|___ANIM.

end

