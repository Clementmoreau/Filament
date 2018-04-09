% Matrix for the 2-filament bundle

function M=matrixNparam_bundle(t,z,N)

global gamma Sp 

zA=z(1:3*N);
zB=z(3*N+1:end);

MA=M_bundle1(t,zA,N);

MA(1,:)=[1,zeros(1,3*N-1)];
MA(2,:)=[zeros(1,N),1,zeros(1,2*N-1)];
MA(3,:)=[zeros(1,2*N),1,zeros(1,N-1)];

MB=M_bundle1(t,zB,N);

MB(1,:)=[1,zeros(1,3*N-1)];
MB(2,:)=[zeros(1,N),1,zeros(1,2*N-1)];
% MB(3,:)=[zeros(1,2*N),1,zeros(1,N-1)];

M=[MA,zeros(3*N);zeros(3*N),MB];


end
