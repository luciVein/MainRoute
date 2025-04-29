clc;
clear;
n = input('Enter the number of varaibles : ');
m = input('Enter the number of constaints : ');
g =input('Enter the number of greater than equation : ');
l =input('Enter the number of less than equation : ');
e =input('Enter the number of euqal equation : ');
A =input('Enter the  matrix : ');
b =input('Enter the b matrix : ');
C = input('Enter the matrix C : ');
con = input('Enter which variable are Integer (if Integer write 1 else 0) :');
con = [con;zeros(e+2*g+l,1)]
con = ones(size(con))-con;
M = 10000;
bvi = zeros(e+g+l,1);
for i =1:g
    t = zeros(m,1);
    t(i,1) = -1;
A = [A t];
C = [C;0];
end
for i =1:g
    t = zeros(m,1);
    t(i,1) = 1;
A = [A t];
C = [C;-M];
bvi(i,1) = n+g+i;
end
for i =1:e
    t = zeros(m,1);
    t(i+g,1) = 1;
A = [A t];
C = [C;-M];
bvi(i+g,1) = n+2*g+i;
end
for i =1:l
    t = zeros(m,1);
    t(i+g+e,1) = 1;
A = [A t];
C = [C;0];
bvi(i+g+e,1) = n+2*g+e+i;
end
A
C
bvi
[m,n] =size(A);
table = zeros(1,n+3);
for i = 4:n+3
    table(1,i) =i-3;
end
c = zeros(m,1);

for i = 1:m
    c(i,1) = C(bvi(i,1),1);
   
end
table = [ table;[bvi c b A]];
table
k = 10000;
while(k>0)
    k = k-1;
    z_c = zeros(1,n);
    for i = 1:n
        z_c(1,i) = (table(:,2)')*table(:,i+3)-C(table(1,i+3),1);
    end
    z_c
    if (any(z_c<0))
        [v,in] = min(z_c);
        pc = in+3;
        y = table(2:m+1,3);
        for i = 1:m
            if(table(i+1,pc)>0)
            y(i,1) = y(i,1)/table(i+1,pc);
            else
                y(i,1) = M;
            end
        end
        [v,in] = min(y);
        if(v==M)
            disp('unbounded solution ');
            return;
        end
        pr = in+1;
        pc
        pr
        ntable = table;
        ntable(pr,1) = pc-3;
        ntable(pr,2) = C(pc-3,1);
        ntable(pr,3:n+3) = table(pr,3:n+3)/table(pr,pc);
        for i = 2:m+1
            if(i~=pr)
                ntable(i,3:n+3) = ntable(i,3:n+3) - table(i,pc)*ntable(pr,3:n+3);
            end
        end
        
        table = ntable;
        table
        continue;
    end
    break;
end
o =3;
while(o)
    y=  table(2:m+1,3);
    y
    z =0;
    for i = 1:m
        
        if(y(i,1) - floor(y(i,1))<1e-5)
            continue;
        end
         if(y(i,1)~=floor(y(i,1))&&con(table(i+1,1),1)==0)
            z = z+1;
        end
    end
        if(z==0)
            disp('Interger solution found ');
            y
            return;
        end
    conv = con;
    for i = 1:size(con,2)
        if con(i,1)==0
            conv(i,1) =2;
        end
    end
    for i =1:m
        conv(table(i+1,1),1) = 0;
    end
    y = mod(y,1);
    y
    [v,in] = max(y);
    in =in+1;
    for i=4:n+3
        if(table(in,i)<0&&conv(i-3,1)==2)
            conv(i-3,1) =3;
        end
    end
    conv
    new = zeros(1,n+3);
    new(1,1)= n+1;
    new(1,2) =0;
    C = [C;0];
    new(1,3) = -1*v;
    for i = 4:n+3
        if(conv(i-3,1)==1)
            new(1,i) = -mod(table(in,i),1);
        end
        if(conv(i-3,1)==2)
            new(1,i) = -table(in,i);
        end
        if(conv(i-3,1)==3)
            new(1,i) = -table(in,i)*v/(1-v);
        end
    end
    table = [table;new];
    t= zeros(m+2,1);
    t(m+2,1) =1;
    t(1,1) = n+1;
    table = [table t];
    m = m+1;
    n = n+1;
    table
    %dual simplex
    k =1000;
    con= [con;1];
   while(k>0)
       k = k-1;
       y= table(2:m+1,3);
       if(any(y<0))
           [v,in] = min(y);
           pr = in+1;
            y=zeros(1,n);
            z_c = zeros(1,n);
         for i = 1:n
             z_c(1,i) = (table(:,2)')*table(:,i+3)-C(table(1,i+3),1);
        
         end
         
         z_c

           for i=4:n+3
               if(table(pr,i)<0)
                  y(1,i-3) = z_c(1,i-3)/table(pr,i);
               else
                   y(1,i-3) = -M;
               end
           end
           [v,in] = max(y);
           pc = in+3;
            pc
            pr
        ntable = table;
        ntable(pr,1) = pc-3;
        ntable(pr,2) = C(pc-3,1);
        ntable(pr,3:n+3) = table(pr,3:n+3)/table(pr,pc);
        for i = 2:m+1
            if(i~=pr)
                ntable(i,3:n+3) = ntable(i,3:n+3) - table(i,pc)*ntable(pr,3:n+3);
            end
        end
        
        table = ntable;
        table
           continue
       end
       break;
   end
   o = o-1;
end
