function [cost, cost1,pg,mcp,util,prodc,pg1,mcp1,util1,prodc1] = optlqg(A,B,C,M1,M2,scale,D,bidtime,var,copies)
cost = 0;cost1=0; options = optimset('Display', 'off');
T= size(C,2);%T= time horizon
M=M1+M2;%M= number of agents
pg = zeros(1,12);mcp=pg;mcp1=mcp;pg1=pg;
M1;%loads
M2;%windmills
% M1 = thermal lads M2 = windfarms
util=0;util1=0;prodc=0;prodc1=0;
%var variance of temperature process
M1=scale*M1; M2=scale*M2; M=scale*M;  MU=zeros(1,M1+M2); SIGMA=blkdiag(eye(M1)*var,eye(M2)*var); X_in=vertcat(ones(M1,1)*(70),ones(M2,1)*100);
n = ones(M,1);K=100;
% X_in=vertcat(ones(M1,1)*70,ones(M2,1)*40);
%dynamics: x(t+1) = A(t) + B(t)U(t) + C(t)
for s=1:copies
    X=X_in;
for t = 1:T
    for user = 1:M1  
        util = util +.5*(X(user)-D(t))^2;
        cost = cost+ .5*(X(user)-D(t))^2; 
    end
   for user = M1+1:M
        cost = cost + .5*X(user)^2;
        prodc = prodc + .5*X(user)^2;
    end
%   U--> 1:T-1 (dim = horem-1) X--> 2:T  (add X(1) also so dim = horem)     
    horem =T-t+1;
    if horem>1
    % define parameters for all users in this horizon
                Aeq_1 = zeros(horem-1,2*horem-1);Aeq_2=Aeq_1;
                for i=1:horem-1
                    Aeq_1(i,i) = -A(1);
                    Aeq_1(i,i+1)=1;
                    Aeq_1(i,i+horem)=-B(1);
                end    
                               
                for i=1:horem-1
                    Aeq_2(i,i) = -A(2);
                    Aeq_2(i,i+1)=1;
                    Aeq_2(i,i+horem)=-B(2);
                end                                                                
             % X(t+1) = AX(t)+BU(t)+C(t)                        
    lambda = ones(horem-1,1);%prices for remaining horizon
    x= zeros(M,2*horem-1);%x(m,1:horem) for state values (x(1) is given), x(m,horem+1:2*horem-1) for inputs   
    %iterations to decide optimal price and input vector
    u = zeros(M,1);
    for priceit = 1:bidtime
        for user = 1:M            
            if user<=M1
%                horem
 %               size(Aeq_1)
  %              size([1,zeros(1,2*horem-2)])
                A1 = vertcat(Aeq_1, [1,zeros(1,2*horem-2)]);
                beq = vertcat(C(1,t:T-1)',X(user));
                H = blkdiag(eye(horem),zeros(horem-1));f = vertcat(-D(t:T)',lambda);
                %lambda = ones(horem,1), purchase cost + desired temp cost
 %               lb = -[Inf;K*ones(horem-1,1);K*ones(horem-1,1)];ub=-lb;
                x(user,:) = quadprog(H,f,[],[],A1,beq,[],[],[],options);%x = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,options)
                u(user)=x(user,horem+1); %stores action at first stage
            else 
                A1 = vertcat(Aeq_2,[1,zeros(1,2*horem-2)]);beq = vertcat(C(2,t:T-1)',X(user));
                H = blkdiag(eye(horem),zeros(horem-1));
                f = vertcat(zeros(horem,1),lambda);
                x(user,:) = quadprog(H,f,[],[],A1,beq,[],[],[],options);%x = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,options)                
 %               x(user,:) = quadprog(H,f,[],[],A1,beq,lb,ub,[],options);%x = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,options)
                %lambda = ones(horem,1) C1 is generation cost
                 u(user)=x(user,horem+1); %stores action at first stage
            end
            
        end 
       lambda = min(subplus(lambda + (1/priceit)*(sum(x(:,horem+1:2*horem-1)))'),K);
 %     lambda = subplus(lambda + (1/priceit)*(sum(x(:,horem+1:2*horem-1)))');
    end

%     1
    % projection
            u = u - ((n'*u)/norm(n)^2)*n;u_im = u;u_im(u_im<0)=0;pg(t) = sum(u_im);
    % state update
  %  u
   % sum(u)
  % lambda'
 
 mcp(t) = lambda(1);
  
  
    for user = 1:M1
        X(user) = A(1)*X(user)+B(1)*u(user) +C(1,t);
    end
    for user = M1+1:M
        X(user) = A(2)*X(user)+B(2)*u(user) +C(2,t);
    end   
    X= X+mvnrnd(MU,SIGMA)';
    end
 %   t
  % X'
end
end
cost=cost/copies;
%implementing myopic policy


for s=1:copies
X=X_in;
    for t = 1:T
    horem =T-t+1;   
    for user = 1:M1
        cost1 = cost1+ .5*(X(user)-D(t))^2;
        util1 = util1 + .5*(X(user)-D(t))^2;
    end
    for user = M1+1:M
        cost1 = cost1+ .5*X(user)^2;
        prodc1 = prodc1+.5*X(user)^2;
    end     
    if horem>1
    lambda = 1;%prices for remaining horizon   
    %iterations to decide optimal price and input vector
    for priceit = 1:bidtime
        u = zeros(M,1);
        for user = 1:M %myopic bid updates
            if user<=M1
%                 priceit
                 u(user) = -(lambda/B(1)+ A(1)*X(user)+C(1,t)-D(t+1))/B(1);  %                (ax+bu+c-d)^2 = 1/2 * (bu)^2 + (bu)(ax+c-d) 
            else
                 u(user) = -(X(user)+lambda);
            end
        %price update
        end
        lambda = min(subplus(lambda + (1/priceit)*sum(u)),K);
    end 
    
%     2
     % projection
             u = u - ((n'*u)/norm(n)^2)*n;u_im = u;u_im(u_im<0)=0;pg1(t) = sum(u_im);mcp1(t) = lambda(1);

        % state update
    for user = 1:M1
        X(user) = A(1)*X(user)+B(1)*u(user) +C(1,t);
    end
    for user = M1+1:M
        X(user) = A(2)*X(user)+B(2)*u(user) +C(2,t);
    end   
    X= X+mvnrnd(MU,SIGMA)';
    end
    end
end
cost1=cost1/copies;
