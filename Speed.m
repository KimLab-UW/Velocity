[filenames_TT_S1,filenames_VT_S1,filenames_TT_S2,filenames_VT_S2,tetrode,Cluster1,Cluster2] = textread('fileinfo.txt','%s%s%s%s%d%d%d');
NN=size(filenames_TT_S1);
N=NN(1);
cellnumber=csvread('fileinfo.csv',0,5);
for k=1:N
    X=filenames_VT_S1{k}
    Y=filenames_TT_S1{k}

    VT = csvread(X,2,0);
    TT = csvread(Y,2,0);
    sc=Cluster1(k,1);
    VT=VT(:,[4 6 7]);
    
    %missing tracking data
    trl = find((VT(:,2)) == 0);
    potx = VT(:,2);
    
    n = 1;
    c = 1; 
    while n <= length(trl)  
        if potx(trl(n)+1) ~= 0         
            up = potx(trl(n) - 1);        
            down = potx(trl(n) + 1);        
            potx(trl(n)) = mean([up down]);        
        elseif potx(trl(n)+1) == 0        
                      
            while potx(trl(n)+c) == 0         
                c = c + 1;        
            end
            up = potx(trl(n) - 1);
            down = potx(trl(n) + c);        
            potx(trl(n):trl(n) + c - 1) = mean([up down]);   
        end
        if c > 1
            n = n + c;        
            c = 1;   
        else
            n = n + 1;
        end
    end
    
    trl = find((VT(:,3)) == 0);
    poty = VT(:,3);
    
    n = 1;
    while n <= length(trl)  
        if poty(trl(n)+1) ~= 0         
            up = poty(trl(n) - 1);        
            down = poty(trl(n) + 1);        
            poty(trl(n)) = mean([up down]);        
        elseif poty(trl(n)+1) == 0        
            c = 1;               
            while poty(trl(n)+c) == 0         
                c = c + 1;        
            end
            up = poty(trl(n) - 1);
            down = poty(trl(n) + c);        
            poty(trl(n):trl(n) + c - 1) = mean([up down]);    
        end
        if c > 1
            n = n + c;        
            c = 1;    
        else
            n = n + 1;
        end
    end
    
    
    VT(:,[2 3])=[potx poty];
    spike=0;
    speedall=[];
    firingrate=[];
    cellTS=TT(find(TT(:,3)==sc),1);
    spike=histc(cellTS,VT(:,1));
      
    sizeVT=size(VT);
    aa=sizeVT(1);
    aa=floor(aa/30);
    for ii=1:aa;
        distpix=sqrt((VT((ii-1)*30+1,2)-VT(ii*30,2))^2 + (VT((ii-1)*30+1,3)-VT(ii*30,3))^2);
        distcm=distpix/2.58;
        time=(VT(ii*30,1)-VT((ii-1)*30+1,1))*0.000001;
        velocity=distcm/time;
        FR=sum(spike((ii-1)*30+1:ii*30))/time;
        
        frate=spike/time;
        speedall=[speedall;velocity];
        firingrate=[firingrate;FR];
    end
    
    temp_s=load('speed_all_S1.txt');
    temp_s=[temp_s,speedall];
    save('speed_all_S1.txt','temp_s','-ascii');
    
    temp_f=load('fr_all_S1.txt');
    temp_f=[temp_f,firingrate];
    save('fr_all_S1.txt','temp_f','-ascii');

end