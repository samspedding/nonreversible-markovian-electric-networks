% Bistochastify a matrix iteratively

function a=bistochastic(a,eps,maxn)
if nargin<2, eps=0.001; end
s1=sum(a,2);
s2=sum(a);
count=0;
while ((abs(max(s1)-1)) | (abs(max(s2)-1)))>eps && count<maxn
  count=count+1;
  a=a./s1;
  s1=sum(a,2);
  s2=sum(a);
  a=a./s2;
end
end
