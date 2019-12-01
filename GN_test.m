function GN_test
U = @(T) 204165.5/(330-2*T) + 10400/(T-20);
%% mínimo U = 1225.2   T= 55.0835
lon=17;
pobla=14;
Timin =40;
Timax =90;
Pmut=0.02;
flag_mut=0;
Iterations=1000;
method=0; % 1 for roulette and  ~1 for Tournament

x=zeros(pobla,lon);
for i=1:pobla
    x(i,:) = randi([0 1],[lon,1])';
end

fprintf('----------------------------------\n')
for itera=1:Iterations
    
    %%%%%%%%% fitness
    f=zeros(pobla,1);
    fenotype=f;
    for i=1:pobla
        DV=bi2de(x(i,:),'left-msb');
        Ti= Timin+ ((Timax-Timin)*DV)/(2^lon-1);
        fenotype(i)=Ti;
        f(i)=U(Ti);
    end
    
    %%%%%%%%%%%%%%%%%%
    [fx_min_curr,ind]=min(f);
    x_min_curr=fenotype(ind);
    if itera==1
        fx_min_best=fx_min_curr;
        x_min_best=x_min_curr;
        fprintf('Minimum = %8.3f  the best fitness fx=%8.6f, x=%8.6f muta= %d \n',1225.2 ,fx_min_best,x_min_best,flag_mut)
    else
        if fx_min_curr<fx_min_best
            fx_min_best=fx_min_curr;
            x_min_best=x_min_curr;
            fprintf('Minimum = %8.3f  the best fitness fx=%8.6f, x=%8.6f muta= %d \n',1225.2 ,fx_min_best,x_min_best,flag_mut)
        else
            fprintf('Minimum = %8.3f  the best fitness fx=%8.6f, x=%8.6f muta= %d \n',1225.2 ,fx_min_best,x_min_best,flag_mut)
        end
    end
    flag_mut=0;
    
    
    
    if method==1
        %%%%%%%%%%%%%%%%%%%%%%%%% Ruleta selección natural
        if min(f)<0
            f=f-min(f);
        end
        
        F=1./(1+f);
        ProF= F/sum(F);
        
        CumP=cumsum(ProF)';
        
        %     %duplicamos la población
        dup = rand([pobla,1])';
        Rep = CumP<dup(:);
        Rep = sum(Rep,2)+1;
        nprogs = x(Rep,:);
        
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%% Selección por Torneo natural
        dup = randi([1 pobla],[pobla,1])';
        % dup=[8,4,2,9,10,7,8,1,4,2];
        
        [~,Tp] = min([f,f(dup)],[],2);
        progs(:,:,1)=x;
        progs(:,:,2)=x(dup,:);
        nprogs=zeros(size(x));
        for i=1:pobla
            nprogs(i,:)=progs(i,:,Tp(i));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Cruce
    pointC=randi([2,lon-1],[1,pobla]);
    ite=1;
    while ite<=pobla
        t=randi([1,pobla],[1,2]);
        if t(1)~=t(2)
            x(ite,1:pointC(ite))=nprogs(t(1),1:pointC(ite));
            x(ite,pointC(ite)+1:lon)=nprogs(t(2),pointC(ite)+1:lon);
            ite=ite+1;
            x(ite,1:pointC(ite-1))=nprogs(t(2),1:pointC(ite-1));
            x(ite,pointC(ite-1)+1:lon)=nprogs(t(1),pointC(ite-1)+1:lon);
            ite=ite+1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Mutación
    mut =find( random('unif',0,1,[1,pobla])<Pmut);
    
    if ~isempty(mut)
        pointm=randi([1,lon],[1,length(mut)]);
        x(mut,pointm)=~x(mut,pointm);
        flag_mut=1;
    end
    
    
end















