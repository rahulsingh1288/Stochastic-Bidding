function[cost, cost1,pg,mcp,util,prodc,pg1,mcp1,util1,prodc1] = tree(A,B,C,P,M1,M2,M3,D,scale,bidtime,c_3,copies,x_in)
% ~ is to be replaced by noise values
%W = [W1 W2; W3 W4]
%copies = 1;
global W
global Buf
cost = 0;cost1=0;pg = zeros(1,2);mcp=pg;pg1=pg;mcp1=mcp;
%P(1,:) contains value of disturbances, P(2,:) contains their probabilities
options = optimset('Display', 'off');
M=M1+M2+M3;%agents number
M;
util=0;util1=0;prodc=0;prodc1=0;  
% Buf is wind farm storage buffer
%M1 loads
%M2 generators
%M3 wind farms
c_1 =.1 ;c_2=.01; % fossil cost = c_1u +c_2u^2+ c_3 (u-x)^2
M1=scale*M1;
M2=scale*M2;
M3=scale*M3;
M=scale*M;
n = ones(M,1);
c_1=-c_1;
Z=vertcat(repmat(A(1,:),M1,1),repmat(A(2,:),M2,1));
A=Z;
Z=vertcat(repmat(B(1,:),M1,1),repmat(B(2,:),M2,1));
B=Z;
Z=vertcat(repmat(C(1,:),M1,1),repmat(C(2,:),M2,1));
C=Z;
X_in=vertcat(ones(M1,1)*x_in(1),x_in(2)*ones(M2,1),x_in(3)*ones(M3,1));  

for s=1:copies 
    X= X_in;    %initial state vector
% P = p11 p12 p21 p22
for time = 1:2
    if time == 1
        % define optimization parameters for thermal loads
                a1 = A(1,1);b1 = B(1,1);c1 = C(1,1); c2 = C(1,2);
                Aeq = zeros(7,10);beq = zeros(7,1);
                Aeq(1,1) = -a1;Aeq(1,2)=1;Aeq(1,8)=-b1;beq(1)=c1;
                Aeq(2,1) = -a1;Aeq(2,3)=1;Aeq(1,8)=-b1;beq(2)=c1;
                Aeq(3,2) = -a1;Aeq(3,4)=1;Aeq(3,9)=-b1;beq(3)=c2; 
                Aeq(4,2) = -a1;Aeq(4,5)=1;Aeq(4,9)=-b1;beq(4)=c2;
                Aeq(5,3) = -a1;Aeq(5,6)=1;Aeq(5,10)=-b1;beq(5)=c2;
                Aeq(6,3) = -a1;Aeq(6,7)=1;Aeq(6,10)=-b1;beq(6)=c2;  
                Aeq(7,1) = 1; 
                %cost 
                H = blkdiag(0, P(1), P(2), P(1)*P(3), P(1)*P(4), P(2)*P(3), P(2)*P(4),zeros(3));%quadratic temp deviation cost
        lambda = ones(3,1);
        u = zeros(M,3);
        for priceit = 1:bidtime
            %generate optimal bid in terms of price vector
            for user = 1:M
            %define cost functions
            if user <= M1 % thermal load
                 beq(7) = X(user) ;
                f = [0 -D(1)*P(1) -D(1)*P(2) -D(2)*P(1)*P(3) -D(2)*P(1)*P(4) -D(2)*P(2)*P(3) -D(2)*P(2)*P(4) lambda(1)  lambda(2)*P(1)  lambda(3)*P(2) ];     
%                f = [0 -D(1)*.5 -D(1)*.5 -D(2)*.25 -D(2)*.25 -D(2)*.25 -D(2)*.25 lambda(1)  lambda(2)*.5  lambda(3)*.5 ];   
                x = quadprog(H,f,[],[],Aeq,beq,[],[],[],options);
 %               priceit
                u(user,1:3) = x(8:10);
            elseif (user > M1) && (user <= M1+M2) %fossil plant
                u(user,:) = fossil(lambda',X(user),c_1,c_2,c_3,P);
            else % wind farm
               u(user,:) = windmill(lambda',X(user),P);
            end
            end
           lambda = subplus(lambda + (1/priceit)*(sum(u))');
        end
       % add fossil generation cost c1*u + c2*u^2 + c3*(u-x)^2 
%      sum(u(:,1))
        mcp(1) = lambda(1);
% projwt u onto sum u_i = 0;
        U = u(:,1);
        U = U - ((n'*U)/norm(n)^2)*n;
        u(:,1) = U;
       for user = M1+1:M1+M2
           cost = cost + c_1*u(user,1)+c_2*u(user,1)^2+c_3*(u(user,1)-X(user))^2;
           prodc = prodc +  c_1*u(user,1)+c_2*u(user,1)^2+c_3*(u(user,1)-X(user))^2;
       end
        U_im = U;U_im(U_im<0)=0;
        pg(1) = pg(1)+sum(U_im);
        
%        % state updates
       n1 = binornd(1,P(1));
       for user = 1:M1+M2
           X(user) = A(user)*X(user) + B(user)*u(user) + C(user,1);
       end      
       for user = M1+M2+1:M
          X(user) = max(min(X(user)+u(user)+W(1,n1+1),Buf),0); 
       end
%        display('opti')
%        u
%        time
%        X'
       % add temperature deviation cost
       for user = 1:M1
           cost = cost + .5*(X(user)-D(1))^2;
           util = util + .5*(X(user)-D(1))^2;
       end
  
    elseif time == 2
        % define optimization parameters for thermal loads   
                beq = zeros(3,1);
                a1 = A(1,2);b1 = B(1,2);c2 = C(1,2);
                Aeq = zeros(3,4);
                Aeq(1,1) = -a1;Aeq(1,2)=1;Aeq(1,4)=-b1;beq(1)=c2;
                Aeq(2,1) = -a1;Aeq(2,3)=1;Aeq(2,4)=-b1;beq(2)=c2;
                Aeq(3,1) = 1;
                H = blkdiag(0,P(3),P(4),0);   
        lambda = 1;
        u = ones(M,1);
        for priceit = 1:bidtime
            %generate optimal bid in terms of price vector
            for user = 1:M
            if user <= M1 % thermal load
                beq(3) = X(user);
                f = [0 -D(2)*P(3) -D(2)*P(4) lambda]; 
                x = quadprog(H,f,[],[],Aeq,beq,[],[],[],options);
                u(user) = x(4);
            elseif (user > M1) && (user <= M1+M2) %fossil plant
                    u(user) = fossil(lambda,X(user),c_1,c_2,c_3,P);
            else % wind farm
                    u(user) = windmill(lambda,X(user),P);
            end
            end
 %           lambda = lambda - (1/priceit)*sum(u);
         lambda = subplus(lambda + (1/priceit)*sum(u));            
        end
        mcp(2) = mcp(2) + lambda;
     %   sum(u)
     % projection
             u = u - ((n'*u)/norm(n)^2)*n; u_im = u;u_im(u_im<0)=0;pg(2) = pg(2)+sum(u_im);
     
     
%      add fossil generation cost c1*u + c2*u^2 + c3*(u-x)^2
       for user = M1+1:M1+M2
           cost = cost + c_1*u(user)+c_2*u(user)^2+c_3*(u(user)-X(user))^2;
           prodc = prodc+  c_1*u(user)+c_2*u(user)^2+c_3*(u(user)-X(user))^2;
       end       
%      state updates
       n2 = binornd(1,.5); 
       for user = 1:M1+M2
           X(user) = A(user)*X(user) + B(user)*u(user) + C(user,2);
       end
       for user = M1+M2+1:M
          X(user) = max(min(X(user)+u(user) + W(2,n2+1),Buf),0); 
       end     
       % add temperature deviation cost
       for user = 1:M1
           cost = cost + .5*(X(user)-D(2))^2;
           util = util + .5*(X(user)-D(2))^2;
       end
%        u
%         time
%        X'
    end    
   
end


% greedy policy
X= X_in;    %initial state vector
for time = 1:2
    if time == 1
                a1 = A(1,1);b1 = B(1,1);c1 = C(1,1); c2 = C(1,2);
                Aeq = zeros(2,3);beq = zeros(2,1);
                Aeq(1,1) = -a1;Aeq(1,2)=1;Aeq(1,3)=-b1;beq(1)=c1;
                Aeq(2,1) = 1;    H = blkdiag(0,1,0);%quadratic temp deviation cost        
        lambda = 1;
        u = zeros(M,1);
        for priceit = 1:bidtime
            %generate optimal bid in terms of price vector
            for user = 1:M
            %define cost functions
            if user <= M1 % thermal load
                % setting parameters
                 beq(2) = X(user) ;
                f = [0 -D(1) lambda]; 
                x = quadprog(H,f,[],[],Aeq,beq,[],[],[],options);
                u(user) = x(3);
            elseif (user > M1) && (user <= M1+M2) %fossil plant
                u(user,:) = fossil(lambda',X(user),c_1,c_2,c_3,P);
            else % wind farm
               u(user,:) = windmill(lambda',X(user),P);
            end
            end

            lambda = subplus(lambda + (1/priceit)*(sum(u))');
        end
        mcp1(1) = mcp1(1) + lambda;
        % projection
        u = u - ((n'*u)/norm(n)^2)*n;u_im = u;u_im(u_im<0)=0;pg1(1) = sum(u_im);

    %     sum(u)
       % add fossil generation cost c1*u + c2*u^2 + c3*(u-x)^2 
       for user = M1+1:M1+M2
           cost1 = cost1 + c_1*u(user,1)+c_2*u(user,1)^2+c_3*(u(user,1)-X(user))^2;
           prodc1 = prodc1 + c_1*u(user,1)+c_2*u(user,1)^2+c_3*(u(user,1)-X(user))^2;
       end
%        % state updates
 
       for user = 1:M1+M2
           X(user) = A(user)*X(user) + B(user)*u(user) + C(user,1);
       end      
       for user = M1+M2+1:M
          X(user) = max(min(X(user)+u(user) + W(1,n1+1),Buf),0); 
       end
    %   X'
       % add temperature deviation cost
       for user = 1:M1
           cost1 = cost1+ .5*(X(user)-D(1))^2;
           util1 = util1+ .5*(X(user)-D(1))^2;
       end
%        u
%        X'
 % second time        
    elseif time == 2
               Aeq = zeros(2,3);beq = zeros(2,1);
                Aeq(1,1) = -a1;Aeq(1,2)=1;Aeq(1,3)=-b1;beq(1)=c2;
                Aeq(2,1) = 1;        
               H = blkdiag(0,1,0);%quadratic temp deviation cost       
        lambda = 1;
        u = ones(M,1);
        for priceit = 1:bidtime
            %generate optimal bid in terms of price vector
            for user = 1:M
            if user <= M1 % thermal load
                beq(2) = X(user) ;
                f = [0 -D(2) lambda]; 
                x = quadprog(H,f,[],[],Aeq,beq,[],[],[],options);
                u(user) = x(3);
            elseif (user > M1) && (user <= M1+M2) %fossil plant
                    u(user) = fossil(lambda,X(user),c_1,c_2,c_3);
            else % wind farm
                u(user) = windmill(lambda,X(user));
            end
            end
            lambda = subplus(lambda + (1/priceit)*sum(u));
        end
         % projection
         mcp1(2)=mcp1(2)+lambda;
         u = u - ((n'*u)/norm(n)^2)*n;u_im = u;u_im(u_im<0)=0;pg1(2) = sum(u_im);

    %     sum(u)
%      add fossil generation cost c1*u + c2*u^2 + c3*(u-x)^2
       for user = M1+1:M1+M2
           cost1 = cost1 + c_1*u(user)+c_2*u(user)^2+c_3*(u(user)-X(user))^2;
           prodc1 = prodc1 + c_1*u(user)+c_2*u(user)^2+c_3*(u(user)-X(user))^2;
       end
       
%      state updates
       for user = 1:M1+M2
           X(user) = A(user)*X(user) + B(user)*u(user) + C(user,2);
       end
       for user = M1+M2+1:M
          X(user) = max(min(X(user)+u(user) + W(2,n2+1),Buf),0); 
       end
    %  X' 
       
       
       % add temperature deviation cost
       for user = 1:M1
           cost1 = cost1 + .5*(X(user)-D(2))^2;
           util1 = util1 +  .5*(X(user)-D(2))^2;
       end
%        u
%        X'
    end    
   
end
end
cost = cost/copies;cost1=cost1/copies;
pg=pg/copies;mcp=mcp/copies;util=util/copies;prodc=prodc/copies;pg1=pg1/copies;mcp1=mcp1/copies;util1=util1/copies;prodc1=prodc1/copies;
