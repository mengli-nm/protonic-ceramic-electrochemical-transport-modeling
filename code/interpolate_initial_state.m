function y=interpolate_initial_state(x)

%%%interpolate_initial_state the first row
sz=size(x);
A1(1:sz(2))=x(1,1:sz(2));
xq1=(1:1:sz(2));
KK1=find(A1>0);
szz1=length(KK1);

for ii=1:1:szz1
    xx1(ii)=KK1(ii);
    yy1(ii)=A1(KK1(ii));
end

vq1 = interp1(xx1,yy1,xq1);
y(1,1:sz(2))=vq1;

%%%interpolate_initial_state the second rows
A2(1:sz(2))=x(2,1:sz(2));
xq2=(1:1:sz(2));
KK2=find(A2>0);
szz2=length(KK2);

for ii=1:1:szz2
    xx2(ii)=KK2(ii);
    yy2(ii)=A2(KK2(ii));
end

vq2 = interp1(xx2,yy2,xq2);
y(2,1:sz(2))=vq2;

%%%interpolate_initial_state the thrid row
A3(1:sz(2))=x(3,1:sz(2));
xq3=(1:1:sz(2));
KK3=find(A3>0);
szz3=length(KK3);

for ii=1:1:szz3
    xx3(ii)=KK3(ii);
    yy3(ii)=A3(KK3(ii));
end

vq3 = interp1(xx3,yy3,xq3);
y(3,1:sz(2))=vq3;

%for jj=2:1:sz(1)
%    y(jj,1:sz(2))=x(jj,1:sz(2));
%end

y(4,1:sz(2))=x(4,1:sz(2));
end