function [soln] = gausselimination( mat1,values,size )
%this function uses the direct method-gaussian elimination to solve a set
%of linear algebraic equations when a matrix is input
soln=zeros(size,1);
beta=values;
for i=1:size
    if i==1	%for the first row of matrix (special case)
        e=mat1(1,1);
    for j=1:size
        mat1(i,j)=mat1(i,j)/e;	%dividing the entire row elements by first element
    end
    beta(1)=beta(1)/e;	
    for k=i+1:size
        c=mat1(k,1);
        for j=1:size
        mat1(k,j)=mat1(k,j)-c*mat1(1,j); %subtracting the elements of the next rows of the first column from the first element
        end
        beta(k)=beta(k)-c*beta(1);
    end
   
    elseif i~=size	%for any other row
        f=mat1(i,i);
        for j=1:size
            mat1(i,j)=mat1(i,j)/f;
        end
        beta(i)=beta(i)/f;

        for k=i+1:size
            d=mat1(k,i);
            for j=1:size
                mat1(k,j)=mat1(k,j)-d*mat1(i,j);
            end
            beta(k)=beta(k)-d*beta(i);
        end

    else
        beta(i)=beta(i)/mat1(i,i);
        mat1(i,i)=1;
    end
end
for m=size:-1:1
    if m==size
        soln(size)=beta(size);
    else
        sum=0.0;
        for g=m+1:size
            sum=sum+mat1(m,g)*soln(g);
        end
        soln(m)=(beta(m,1)-sum)/mat1(m,m);	%calculating the solution by backward sweeping
    end
    
 end