%计算模拟ts，xs即为模拟ts
function [xs,nn]=simulAGint(a,G,C,long_total,w,TRsec,val)%C是SC %val是重采样的时间点

    
    %RECORDAR QUE ESTE RESOLVEDOR PERMITE TENER UN VECTOR DE G, Y NO ES
    %IGUALAL OTRO YA QUE SACA G Y LO INTRODUCE DIRECTO EN LA EC.
    %POR AHORA EL G ES EL G INTERNO
    dt = 0.6;
    sig = 0.04;%.04; %was 0.04
    dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step
    
    nNodes=length(a);
    
    a=repmat(a,1,2);
    
    wC = C ;%aca es donde saco el G afuera de como estaba antes %这是我把G拿出来的地方，就像之前那样
    
    wintC = wC.*repmat(G',nNodes,1);%SC*一个90*90的值都是0.5的矩阵，意思是SC与G结合
    
    G=repmat(G,1,2);


    A = 0;%扰动幅值
    region = 15;   %扰动区域 
    
     

    xs=zeros(nNodes,long_total);
    z = 0.1*ones(nNodes,2); % --> x = z(:,1), y = z(:,2)
    zz = z(:,end:-1:1);
    
    sumC = repmat(wC*G(:,1),1,2);

    %comienzo a simular con el a otorgado (RECORDAR QUE ACA a YA NO ES
    %HOMOGENEO Y TAMPOCO G



        nn=1;
        
        for t=1:dt:1000 %Ignacio estuvo aqu�, y bajo el tiempo de termalizacion.
            %suma =  wC*(z.*G) - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            
           zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            
           z = z + dt*(a.*z + zz.*w - z.*(z.*z+zz.*zz) + wintC*(z) - sumC.*z)  + dsig*randn(nNodes,2);
           
        end
        
        t=1:dt:long_total*TRsec;
      
       
        for i=1:length(t)
            %suma = wC*(z.*G) - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
           
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*w - z.*(z.*z+zz.*zz) + wintC*(z) - sumC.*z)  + dsig*randn(nNodes,2); 

            z(region,1)=z(region,1) + A*cos(w(region,2)*(i-1)*dt);
            z(region+1,1)=z(region+1,1) + A*cos(w(region+1,2)*(i-1)*dt);
            z(region,2)=z(region,2) + A*sin(w(region,2)*(i-1)*dt);
            z(region+1,2)=z(region+1,2) + A*sin(w(region+1,2)*(i-1)*dt);
            
            if t(i)==val(nn)
                
                %RECORDAR: por alguna razon se les ocurrio poner las tseries en
                %columnas. Yo lo corrijo para no perder sanidad
                xs(:,nn)=z(:,1);
                nn=nn+1;
                
            end
        
        end
        
        