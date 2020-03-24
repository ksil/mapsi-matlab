function val = cylinderformfactor(R,L,in,out)

phi = in(:,1);
theta = in(:,2);
Qx = out(:,1);
Qy = out(:,2);
Qz = out(:,3);

ct = cos( theta );
st = sin( theta );
cp = cos( phi );
sp = sin( phi );

tmp1 = R.*sqrt((Qy.*cp - Qx.*sp).^2+(ct.*(Qx.*cp+Qy.*sp) - Qz.*st).^2);
tmp2 = (1/2).*L.*(Qz.*ct+(Qx.*cp+Qy.*sp).*st);

val = 2.*jinc(tmp1).*mysinc(tmp2);
val = val.*val;

end

function val = jinc(x)

val = besselj(1,x) ./ x;
val(x == 0.0) = 0.5;

end

function val = mysinc(x)

val = sin(x) ./ x;
val(x == 0.0) = 0.0;

end