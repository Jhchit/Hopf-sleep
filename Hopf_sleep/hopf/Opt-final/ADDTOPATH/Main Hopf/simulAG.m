
function [xs]=simulAG(a,G,C,long_total,w,TRsec,val)

    
    %RECORDAR QUE ESTE RESOLVEDOR PERMITE TENER UN VECTOR DE G, Y NO ES
    %IGUAL AL OTRO YA QUE SACA G Y LO INTRODUCE DIRECTO EN LA EC.
    %POR AHORA EL G ES EL G ASIMETRICO
    %“%记住，这个求解器允许有一个G向量，且与另一个不同，因为它将G提取并直接引入方程中。%目前G是非对称的G。”
    dt = 0.1;
    sig = 0.04; %was 0.04
    dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step
    
    nNodes=length(a);
    
    a=repmat(a,1,2);
    
    G=repmat(G,1,2);
        
    wC = C;%aca es donde saco el G afuera de como estaba antes 

    
    xs=zeros(nNodes,long_total);
    z = 0.1*ones(nNodes,2); % --> x = z(:,1), y = z(:,2)

    sumC = repmat(sum(wC,2),1,2);

    %comienzo a simular con el a otorgado (RECORDAR QUE ACA a YA NO ES
    %HOMOGENEO Y TAMPOCO G
     %"%开始用所给予的a进行模拟（记住，这里的a不再是均匀（这个可能是指a的分组）的，G也不是）。"



        nn=1;
        for t=1:dt:1000 %JVS is it really necessary to swing in for 3000secs?
            %suma =  wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*w - z.*(z.*z+zz.*zz) + G.* (wC*z - sumC.*z)) + dsig*randn(nNodes,2);
        end
        
        t=1:dt:long_total*TRsec;
        

        for i=1:length(t) %JVS: was 15000, now faster
            %suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*w - z.*(z.*z+zz.*zz) + G.* (wC*z - sumC.*z)) + dsig*randn(nNodes,2);

            if t(i)==val(nn)
                
                %RECORDAR: por alguna razon se les ocurrio poner las tseries en
                %columnas. Yo lo corrijo para no perder sanidad
                xs(:,nn)=z(:,1);
                nn=nn+1;
            end
        end
        
end

