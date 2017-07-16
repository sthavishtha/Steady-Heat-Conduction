function [ t_old,t_new,diff] = line_by_liney( t_old,t_new,z1,z2,z3,m,n,value,mat,diff )
%this function implements the line by line gauss siedel method along the
%y-direction
mat1=zeros(n,n);
value1=zeros(1,n);
for j=z1:z2:z3
    if j==1
        for i=1:m
            if i==1
                value1(i)= (mat(i,2)*t_old(2,1)-value(1));
                mat1(i,i+1)=mat(j,i+n);
                mat1(j,j)=mat(j,j);
            elseif i==m
                value1(i)=(mat(j+(i-1)*n,j+(i-1)*n+1)*t_old(j+(i-1)*n+1,1)-value(j+(i-1)*n));
                mat1(i,i-1)=mat(j+(i-1)*n,j+(i-2)*n);
                mat1(i,i)=mat(j+(i-1)*n,j+(i-1)*n);
            else
                value1(i)=(mat(j+(i-1)*n,j+(i-1)*n+1)*t_old(j+(i-1)*n+1,1)-value(j+(i-1)*n));
                mat1(i,i+1)=mat(j+(i-1)*n,j+(i)*n);
                mat1(i,i)=mat(j+(i-1)*n,j+(i-1)*n); 
                mat1(i,i-1)=mat(j+(i-1)*n,j+(i-2)*n);
            end
        end
            t_new(1:n:(m-1)*n+1,1)=tdma(mat1,-value1,m);
            for k=1:n:(m-1)*n+1
                diff1=abs(t_new(k,1)-t_old(k,1));
                diff=max(diff,diff1);
            end
            t_old(1:n:(m-1)*n+1,1)=t_new(1:n:(m-1)*n+1,1);
    elseif j==m
        for i=1:m
            if i==1
                value1(i)=(mat(j+(i-1)*n,j+(i-1)*n-1)*t_old(j+(i-1)*n-1,1)-value(j+(i-1)*n));
                mat1(i,i)=mat(j+(i-1)*n,j+(i-1)*n);    
                mat1(i,i+1)=mat(j+(i-1)*n,j+i*n);
            elseif i==m
                value1(i)=(mat(j+(i-1)*n,j+(i-1)*n-1)*t_old(j+(i-1)*n-1,1)-value(j+(i-1)*n));
                mat1(i,i-1)=mat(j+(i-1)*n,j+(i-2)*n);
                mat1(i,i)=mat(j+(i-1)*n,j+(i-1)*n);
            else
                value1(i)=(mat(j+(i-1)*n,j+(i-1)*n-1)*t_old(j+(i-1)*n-1,1)-value(j+(i-1)*n));
                mat1(i,i-1)=mat(j+(i-1)*n,j+(i-2)*n);
                mat1(i,i)=mat(j+(i-1)*n,j+(i-1)*n);    
                mat1(i,i+1)=mat(j+(i-1)*n,j+i*n);
            end
        end
            t_new(m:m:m*n,1)=tdma(mat1,-value1,m);
            for k=m:m:m*n
                diff1=abs(t_new(k,1)-t_old(k,1));
                diff=max(diff,diff1);
            end
            t_old(m:m:m*n,1)=t_new(m:m:m*n,1);
    else
         for i=1:m
            if i==1
                value1(i)=(mat(j+(i-1)*n,j+(i-1)*n-1)*t_old(j+(i-1)*n-1,1)-value(j+(i-1)*n)+mat(j+(i-1)*n,j+(i-1)*n+1)*t_old(j+(i-1)*n+1,1));
                mat1(i,i)=mat(j+(i-1)*n,j+(i-1)*n);    
                mat1(i,i+1)=mat(j+(i-1)*n,j+i*n);
            elseif i==m
                value1(i)=(mat(j+(i-1)*n,j+(i-1)*n-1)*t_old(j+(i-1)*n-1,1)-value(j+(i-1)*n)+mat(j+(i-1)*n,j+(i-1)*n+1)*t_old(j+(i-1)*n+1,1));
                mat1(i,i-1)=mat(j+(i-1)*n,j+(i-2)*n);
                mat1(i,i)=mat(j+(i-1)*n,j+(i-1)*n);
            else
                value1(i)=(mat(j+(i-1)*n,j+(i-1)*n-1)*t_old(j+(i-1)*n-1,1)-value(j+(i-1)*n)+mat(j+(i-1)*n,j+(i-1)*n+1)*t_old(j+(i-1)*n+1,1));
                mat1(i,i-1)=mat(j+(i-1)*n,j+(i-2)*n);
                mat1(i,i)=mat(j+(i-1)*n,j+(i-1)*n);    
                mat1(i,i+1)=mat(j+(i-1)*n,j+i*n);
            end
         end       
    end
            t_new(j:m:j+(m-1)*m,1)=tdma(mat1,-value1,m);
            for k=j:m:j+(m-1)*m
                diff1=abs(t_new(k,1)-t_old(k,1));
                diff=max(diff,diff1);
            end
            t_old(j:m:j+(m-1)*m,1)=t_new(j:m:j+(m-1)*m,1);
end
end

