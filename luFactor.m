function [L,U,P] = luFactor(A)
%luFactor: determines the LU factorization of a square matrix using;
%partial pivoting;
% Inputs/Outputs;
%A = coefficient matrix;
%L= lower triangle matrix;
%U= upper triangle matrix;
%P= the pivot matrix;
format long
if nargin<1, error('coefficient matrix is required')
end
[n,m]=size(A); %n is # of rows in A, m is # of columns in A
if n~=m, error('matrix A must be a square matrix')
end
n=length(A);
P=eye(n)
L=P

for i=[2:n]
x=A(:,1)
[val,loc]=max(abs(x)) %determine which row contain the pivot element 
d=loc                 %d=the row that will be pivoted if necessary
r=A(1,:)  %first row of OG matrix A
s=P(1,:)  %fist row of identity matrix
A(1,:)=A(d,:); %swaps R1 of A with row with row with the pivot element 
A(d,:)=r      % replaces the row with the pivot element with OG row 1 from matrix A
P=P
P(1,:)=P(d,:); %swaps R1 of P with the row with the pivot element
P(d,:)=s      %replaces the row with the pivot element with OG row 1 from matrix P
end
for i=2:n       
    Li1=A(i,1)/A(1,1)    %ratio needed to eleminate the first element of r2
    L(i,1)=Li1 %value used to eleminate the first element of r2, first value for L matrix
    A(i,:)=A(i,:)-Li1*A(1,:) %subtract modified row 1 from row 2
    
end
for i=3:n
    y=A(:,2) %selects column 2
    [val,loc]=max(abs(y)) 
    e=loc
    t=A(2,:) %second row of NEW matrix A
    v=P(2,:) %second row of NEW matrix P
    A(2,:)=A(e,:); %swaps R2 of NEW A with row with the pivot element
    A(e,:)=t %replaces the row that had the pivot element with second row of OLD matrix A
    P=P
    P(2,:)=P(e,:); %swaps R2 of NEW P with the row with the pivot element
    P(e,:)=v      %replaces the row with the pivot element with row 2 from old matrix P
    
end
for i=3:n
    Li3=A(i,2)/A(2,2)
    L(i,2)=Li3
    A(i,:)=A(i,:)-Li3*A(2,:)
end
for i=4:n
    z=A(:,3)%selects column 3
    [val,loc]=max(abs(z))
    z2=[0; 0; 1; 1;].*z %corrects the location of the highest value in the 3rd column
    [val,loc]=max(abs(z2))
    f=loc
    w=A(3,:) %third row of NEW matrix A
    g=P(3,:) %third row of NEW matrix P
    A(3,:)=A(f,:); %swaps R3 of NEW A matrix with row with the pivot element
    A(f,:)=w %replaces the row that had the pivot element with the third row of OLD matrix A
    P=P
    P(3,:)=P(f,:); %swaps R3 of NEW P with the row with the pivot element
    P(f,:)=g %replaces the row with the pivot element 
end
for i=4:n
    Li4=A(i,3)/A(3,3)
    L(i,3)=Li4
    A(i,:)=A(i,:)-Li4*A(3,:)
end
U=A
L=L
    
