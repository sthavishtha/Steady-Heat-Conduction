mat=zeros(2000,200);
value=zeros(1,2000);
alpha=1.25;
disp('steady state heat conduction problem');
m=input('enter the number of nodes in the x direction');
n=input('enter the number of nodes in y direction');
t1=input('enter the temperature of left wall');
t2=input('enter the temperature of right wall');
t3=input('enter the temperature of bottom wall');
t4=input('enter the temperature of top wall');
if m==''
    m=20;
elseif n==''
    n=20;
elseif t1==''
    t1=300;
elseif t2==''
    t2=100;
elseif t3==''
    t3=400;
elseif t4==''
    t4=200;
end
t=zeros(m,n);
t_old=zeros(m*n,1);
t_new_d=zeros(m*n,1);
t_new=zeros(m*n,1);
deltax=1/m;
deltay=1/n;
w1=1/(deltax*deltax);
w2=1/(deltay*deltay);

for i=1:m
    for j=1:n
        if i==1 
            if j==1 
                ap=-3*(w1+w2);
                ae=w1;
                an=w2;
                s=2*t1*w1+2*t3*w2;
                value(1)=-s;
                for k=1:m
                    for l=1:n
                        if k==1 && l==1
                            mat(1,1)=ap;
                        elseif k==2 && l==1
                            mat(1,n+1)=an;
                        elseif k==1 && l==2
                            mat(1,2)=ae;
                        end
                    end
                end
            elseif j==n
                ap=-3*(w1+w2);
                aw=w1;
                an=w2;
                s=2*t2*w1+2*t3*w2;
                value(j)=-s;
                for k=1:m
                    for l=1:n
                        if k==1 && l==n
                            mat(n,n)=ap;
                        elseif k==1 && l==n-1
                            mat(n,n-1)=aw;
                        elseif k==2 && l==n
                            mat(n,2*n)=an;
                        end
                    end
                end
            else
                ap=-(2*w1+3*w2);
                s=2*t3*w2;
                ae=w1;
                aw=w1;
                an=w2;
                value(j)=-s;
                for k=1:m
                    for l=1:n


                                if k==1 && l==j
                                    mat(j+(i-1)*n,l+(k-1)*n)=ap;  
                                elseif k==1 && l==j-1
                                    mat(j+(i-1)*n,l+(k-1)*n)=aw;
                                elseif k==1 && l==j+1
                                    mat(j+(i-1)*n,l+(k-1)*n)=ae;
                                elseif k==2 && l==j
                                    mat(j+(i-1)*n,l+(k-1)*n)=an;
                                end

                    end
                end
            end
        elseif i==m 
            if j==1
                ap=-3*(w2+w1);
                as=w2;
                ae=w1;
                s=2*(t1*w1+t4*w2);
                value(j+(i-1)*n)=-s;
                for k=1:m
                    for l=1:n
                        if k==m && l==1
                            mat(j+(i-1)*n,l+(k-1)*n)=ap;
                        elseif k==m-1 && l==1
                            mat(j+(i-1)*n,l+(k-1)*n)=as;
                        elseif k==m && l==2
                            mat(j+(i-1)*n,l+(k-1)*n)=ae;
                        end
                    end
                end
            elseif j==n
                ap=-3*(w2+w1);
                as=w2;
                aw=w1;
                s=2*(t2*w1+t4*w2);
                value(j+(i-1)*n)=-s;
                for k=1:m
                    for l=1:n
                        if k==m && l==n
                            mat(j+(i-1)*n,l+(k-1)*n)=ap;
                        elseif k==m && l==n-1
                            mat(j+(i-1)*n,l+(k-1)*n)=aw;
                        elseif k==m-1 && l==n
                             mat(j+(i-1)*n,l+(k-1)*n)=as;
                        end
                    end
                end
            else
                ap=-(2*w1+3*w2);
                s=2*w2*t4;
                ae=w1;
                aw=w2;
                value(j+(i-1)*n)=-s;
                for k=1:m
                    for l=1:n

                            if k==m && l==j
                               mat(j+(i-1)*n,l+(k-1)*n)=ap;  
                            elseif k==m-1 && l==j 
                                 mat(j+(i-1)*n,l+(k-1)*n)=as;
                            elseif k==m && l==j-1 
                                 mat(j+(i-1)*n,l+(k-1)*n)=aw;
                            elseif k==m && l==j+1 
                                 mat(j+(i-1)*n,l+(k-1)*n)=ae;
                            end

                       
                    end
                end
            end
        elseif j==1
            if i~=1 && i~=m
                ap=-(3*w1+2*w2);
                s=2*t1*w1;
                ae=w1;
                an=w2;
                as=w2;
                value(j+(i-1)*n)=-s;
                for k=1:m
                    for l=1:m
                        if k==i && l==j
                            mat(j+(i-1)*n,l+(k-1)*n)=ap;
                        elseif k==i-1 && l==j
                             mat(j+(i-1)*n,l+(k-1)*n)=as;
                        elseif k==i+1 && l==j
                             mat(j+(i-1)*n,l+(k-1)*n)=an;      
                        elseif k==i && l==j+1
                            mat(j+(i-1)*n,l+(k-1)*n)=ae;
                        end
                    end
                end
            end
        elseif j==n
             if i~=1 && i~=m
                 ap=-(3*w1+2*w2);
                 s=2*t2*w1;
                 aw=w1;
                 an=w2;
                 as=w2;
                 value(j+(i-1)*n)=-s;
                 for k=1:m
                     for l=1:n
                         if k==i && l==j
                             mat(j+(i-1)*n,l+(k-1)*n)=ap;
                         elseif k==i-1 && l==j
                             mat(j+(i-1)*n,l+(k-1)*n)=as;
                          elseif k==i+1 && l==j
                             mat(j+(i-1)*n,l+(k-1)*n)=an;
                         elseif k==i && l==j-1
                             mat(j+(i-1)*n,l+(k-1)*n)=aw;
                         end
                     end
                 end
             end
        else
            ap=-2*(w1+w2);
            ae=w1;
            aw=w1;
            an=w2;
            as=w2;
            s=0.0;
            value(j+(i-1)*n)=-s;
            for k=1:m
                for l=1:n
                    if k==i && l==j
                         mat(j+(i-1)*n,l+(k-1)*n)=ap;
                    elseif k==i-1 && l==j
                             mat(j+(i-1)*n,l+(k-1)*n)=as;
                    elseif k==i+1 && l==j
                             mat(j+(i-1)*n,l+(k-1)*n)=an;
                    elseif k==i && l==j-1
                             mat(j+(i-1)*n,l+(k-1)*n)=aw;
                    elseif k==i && l==j+1
                        mat(j+(i-1)*n,l+(k-1)*n)=ae;
                    end
                end
            end
        end
    end
end

iter=0;
for i=1:m
    for j=1:n
     t_old(j+(i-1)*n,1)=min(t1,t3); %temperature at all nodal points is initialised 
    end
end
diff=40;
n_d=1;
rat=(1-alpha)/alpha;
while diff>1e-02
    iter=n_d;
for i=1:m*n     %moving through all nodal points and finding temperature
 sum=0.0;

    if i==1
        iter=1*n;
 
        for k=1:m*n
            if k~=i
        sum=sum+mat(i,k)*t_old(k);
            end
        end
        t_new_d(i)=(value(i)-sum)/mat(1,1);
        ra=mat(1,1)/alpha;
        t_new(i)=((t_new_d(i)*mat(1,1))+(rat*mat(1,1)*t_old(1)))/ra;
        
    diff=abs(t_old(1)-t_new(1));
    t_old(1)=t_new(1);
    else

         for k=1:m*n
             if k~=i
                 sum=sum+mat(i,k)*t_old(k);
             end
         end
        t_new_d(i)=(value(i)-1*sum)/mat(i,i);
        ra=mat(i,i)/alpha;
        t_new(i)=((t_new_d(i)*mat(i,i))+(rat*mat(i,i)*t_old(i)))/ra;
     end
     diff1=abs(t_old(i)-t_new(i));
     t_old(i)=t_new(i);
     diff=max(diff,diff1);   
end
    n_d=n_d+1;
end
for s=1:m
    for q=1:n
        t(s,q)=t_new(q+(s-1)*n);
    end
end
[x,y]=meshgrid(1:m,1:n);
contourf(x,y,t);
colorbar;