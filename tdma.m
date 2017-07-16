function [ soln ] = tdma( mat1,values,size )
%this function exhibits the tri-diagonal matrix algorithm, also called as
%thomas algorithm
%this is a direct method of solving a system of linear equations and the
%solution matrix obtained by this function has been validated with some
%standard cases given in Veertseg book
gamma=zeros(size,1);
beta=zeros(size,1);
soln=zeros(size,1);
for i=1:size
    if i==1
        gamma(1,1)=mat1(1,2)/mat1(1,1);
        beta(1,1)=values(1,1)/mat1(1,1);
    else
        if i~=size %gamma does not exist for the last row of elements
        dr=mat1(i,i)-mat1(i,i-1)*gamma(i-1,1);
        gamma(i,1)=mat1(i,i+1)/dr;
        nr=values(i)-mat1(i,i-1)*beta(i-1,1);
        beta(i,1)=nr/dr;
        else
         dr=mat1(i,i)-mat1(i,i-1)*gamma(i-1,1);
         nr=values(i)-mat1(i,i-1)*beta(i-1,1);
        beta(i,1)=nr/dr;
        end
    end
end
   for i=size:-1:1  %calculating the solution matrix by bcakward sweeping 
       if i==size
           soln(size,1)=beta(size,1);
       else
           soln(i,1)=beta(i,1)-gamma(i,1)*soln(i+1,1);
       end
   end
  
end

