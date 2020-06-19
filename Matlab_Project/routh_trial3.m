%% routh hurwitz criteria
clear all
close all
clc
%% firstly it is required to get first two row of routh matrix

e=input('1 3 13 23 12 0 ');
disp('-------------------------------------------------------------------------')
L=length(e);
m=mod(L,2);
if m==0
    a=rand(1,(L/2));  %return matrix
    b=rand(1,(L/2));
    for i=1:(L/2)
        a(i)=e((2*i)-1);
        b(i)=e(2*i);
    end
else
    e1=[e 0];
    a=rand(1,((L+1)/2));
    b=[rand(1,((L-1)/2)),0];
    for i=1:((L+1)/2)
        a(i)=e1((2*i)-1);
        b(i)=e1(2*i);
    end
end

%% now we generate the remaining rows of routh matrix

L1=length(a);
c=zeros(L,L1);
c(1,:)=a;
c(2,:)=b;
for m=3:L
    for n=1:L1-1
        c(m,n)=-(1/c(m-1,1))*det([c((m-2),1) c((m-2),(n+1));c((m-1),1) c((m-1),(n+1))]);
    end
    % Auxillary Form
    if c(m,:) == 0
        ax = zeros(1,L-m+2);
        d = 2;
        for p=1 : L-m+2
            if mod(p,2)
                if p==1
                    ax(p) = c(m-1,p);
                else
                    ax(p) = c(m-1,d);
                    d = d + 1;
                end
            else
                ax(p) = 0;
            end
        end
        z = polyder(ax);
    d = 2;
    for q=1 : length(z)
        if mod(q,2)
            if q==1
                c(m,q) = z(q);
            else
                c(m,d) = z(q);
                d = d+1;
            end
        end
    end
    % Epsilon Form
    else
        if c(m,1) == 0
            c(m,1) = 0.01
        end
    end
           
end
disp('the routh matrix:')
disp(c)
%% now we check the stablity of system
if c(:,1)>0
    disp('System is Stable')
else
    disp('System is Unstable');
end
