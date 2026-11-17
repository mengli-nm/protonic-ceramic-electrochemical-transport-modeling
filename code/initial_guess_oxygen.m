function x0=initial_guess_oxygen(Ph2o,Po2,Cd)
xL=0.0+1e-3;
xR=1.0+1e-3;
x0=0.0;
IterR=1000;

for kk=1:1:10
    xx=(xL:(xR-xL)/IterR:xR);
    for ii=1:1:length(xx)
        yy(ii)=charge_balance_oxygen(xx(ii),Ph2o,Po2,Cd);
    end
    cs = yy.*circshift(yy,-1,2); 
%plot(xx,yy)  
%    Ind0=find(cs==0);
%    if (~isempty(Ind0))
%      if (yy(Ind0(1))==0)
%          x0=xx(Ind0(1));
%      elseif (yy(Ind0(2))==0)
%          x0=xx(Ind0(2));
%      end
%        break;
%    end
         
   Ind=find(cs<=0);
   if (~isempty(Ind))
       if (yy(Ind(1))==0)
         x0=xx(Ind(1));
         break;
       elseif (yy(Ind(2))==0)
         x0=xx(Ind(2));
         break;
       else
         xL=xx(Ind(1));
         xR=xx(Ind(2));
       end
   else
      fprintf('game over \n');
      break;
   end
end
if (x0==0.0)
   x0=xL;
end
end