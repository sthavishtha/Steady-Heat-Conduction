function [ t_old,t_new,diff ] = line_by_linex( t_old,t_new,z1,z2,z3,m,n,value,mat,diff )
%this function implements the line by line gauss siedel method along the
%x-direction
mat1=zeros(n,n);
value1=zeros(1,n);
for i=z1:z2:z3
  if i==1
    for j=1:n

            if j==1
           value1(j)= (mat(i,n+j)*t_old(j+(i)*n,1)-value(j+(i-1)*n));
           mat1(j,j+1)=mat(i,j+1);
           mat1(j,j)=mat(i,i);
            elseif j==n
                value1(j)=(mat(j,n+j)*t_old(j+(i)*n,1)-value(j+(i-1)*n));
                mat1(j,j-1)=mat(j+(i-1)*n,j-1);
                mat1(j+(i-1)*n,j)=mat(j+(i-1)*n,j);
            else
                value1(j)=(mat(j,n+j)*t_old(j+i*n,1)-value(j+(i-1)*n));
                mat1(j,j-1)=mat(j,j-1);
                mat1(j,j+1)=mat(j,j+1);
                 mat1(j+(i-1)*n,j)=mat(j+(i-1)*n,j);
            end
            t_new(1:n,1)=tdma(mat1,-value1,n);
            diff=abs(t_new(1,1)-t_old(1,1));
            for k=2:n
                diff1=abs(t_new(k,1)-t_old(k,1));
                diff=max(diff,diff1);
            end
            t_old(1:n,1)=t_new(1:n,1);
            
    end
  elseif i==m
      for j=1:n
                if j==1
                    value1(j)=(mat(j+(i-1)*n,j+(i-2)*n)*t_old(j+(i-2)*n,1)-value(j+(i-1)*n));
                    mat1(j,j+1)=mat(j+(i-1)*n,j+(i-1)*n+1);
                    mat1(j,j)=mat(j+(i-1)*n,j+(i-1)*n);
                
                elseif j==n
                    value1(j)=(mat(j+(i-1)*n,j+(i-2)*n)*t_old(j+(i-2)*n,1)-value(j+(i-1)*n));
                    mat1(j,j-1)=mat(j+(i-1)*n,j+(i-1)*n-1);
                    mat1(j,j)=mat(j+(i-1)*n,j+(i-1)*n);
                else
                    value1(j)=(mat(j+(i-1)*n,j+(i-2)*n)*t_old(j+(i-2)*n,1)-value(j+(i-1)*n));
                    mat1(j,j-1)=mat(j+(i-1)*n,j+(i-1)*n-1);
                    mat1(j,j+1)=mat(j+(i-1)*n,j+(i-1)*n+1);
                    mat1(j,j)=mat(j+(i-1)*n,j+(i-1)*n);
                end
              t_new((m-1)*n+1:m*n,1)=tdma(mat1,-value1,n);
              for k=(m-1)*n+1:m*n
                diff1=abs(t_new(k,1)-t_old(k,1));
                diff=max(diff,diff1);
              end
              t_old((m-1)*n+1:m*n,1)= t_new((m-1)*n+1:m*n,1);
      end
  else
      for j=1:m
          if j==1
              value1(j)=(mat(j+(i-1)*n,j+(i-2)*n)*t_old(j+(i-2)*n,1)-value(j+(i-1)*n)+mat(j+(i-1)*n,j+i*n)*t_old(j+i*n,1));
              mat1(j,j+1)=mat(j+(i-1)*n,j+(i-1)*n+1);
              mat1(j,j)=mat(j+(i-1)*n,j+(i-1)*n);    
          elseif j==m
              value1(j)=(mat(j+(i-1)*n,j+(i-2)*n)*t_old(j+(i-2)*n,1)-value(j+(i-1)*n)+mat(j+(i-1)*n,j+i*n)*t_old(j+i*n,1));
              mat1(j,j-1)=mat(j+(i-1)*n,j+(i-1)*n-1);
              mat1(j,j)=mat(j+(i-1)*n,j+(i-1)*n); 
          else
              value1(j)=(mat(j+(i-1)*n,j+(i-2)*n)*t_old(j+(i-2)*n,1)-value(j+(i-1)*n)+mat(j+(i-1)*n,j+i*n)*t_old(j+i*n,1));
              mat1(j,j-1)=mat(j+(i-1)*n,j+(i-1)*n-1);
              mat1(j,j+1)=mat(j+(i-1)*n,j+(i-1)*n+1);
              mat1(j,j)=mat(j+(i-1)*n,j+(i-1)*n);           
          end
            t_new((i-1)*n+1:(i)*n,1)=tdma(mat1,-value1,n);
            for k=(i-1)*n+1:(i)*n
                diff1=abs(t_new(k,1)-t_old(k,1));
                diff=max(diff,diff1);
            end    
            t_old((i-1)*n+1:(i)*n,1)=t_new((i-1)*n+1:(i)*n,1);
      end
  end
end
end

