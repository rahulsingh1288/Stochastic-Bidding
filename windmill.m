function[bid] = windmill(lambda,level,P)
options = optimset('Display', 'off');
global Buf
t = size(lambda,2);
%p1=1;p2=0;
if t==3
p1=P(1);p2=P(2);
% x = (x1,x2,x3,u1,u2,u3) not confuse with x in tree.m
fun = @(x)lambda(1)*x(4)+p1*lambda(2)*x(5)+p2*lambda(3)*x(6);
%fun = @(x)-lambda(1)*x(4)-p1*lambda(2)*x(5)-p2*lambda(3)*x(6);
% supply means -ve; so  -level<=u
% withdraw means +ve, so u<= B-level
A = zeros(3+3,6);
A(1,1) = -1;A(1,4)=-1;A(2,2)=-1;A(2,5)=-1;A(3,3)=-1;A(3,6)=-1;% -x<=u
A(4,1) = 1;A(4,4) =1;A(5,2)=1;A(5,5)=1;A(6,3)=1;A(6,6)=1;  % u<=B-x                             
b=zeros(3+3,1);b(4) = Buf;b(5)=Buf;b(6)=Buf;
Aeq = [1 0 0 0 0 0]; beq = level;x0=ones(1,6);
%display(A)
x = fmincon(fun,x0,A,b,Aeq,beq,[],[],@storagecon2,options);

bid = x(4:6);
elseif t==1
 %  x = (x1,x2,x3,u1)
% supply means -ve; so  -level<=u
% withdraw means +ve, so u<= B-level

 if lambda >0
     bid = -level;
 else
     bid = (Buf-level);         
 end
%  if lambda >0
%     bid = level;
%  else
%      bid = -(Buf-level);
%  end
 
end
end

function[c,ceq] = storagecon2(x)
global W
global Buf
n1=W(1,1);n2=W(1,2);
ceq(1) = max(min(x(1)+x(4)+n1,Buf),0)-x(2); % x2 = x1+u
ceq(2) = max(min(x(1)+x(4)+n2,Buf),0)-x(3);
c = [];
end
