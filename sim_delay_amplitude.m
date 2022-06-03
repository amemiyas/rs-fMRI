function sim_delay_amplitude
x1=textread('C:\Users\Precision\Desktop\t1.txt','%f');
nframes=1200;
TR=0.72;
y1 = wgn(nframes*10*100,1,20);
z1 = wgn(nframes*10*100,1,20);
R=zeros(60,60,100,3);
mR=zeros(60,60,3);

for k=1:3 % noise level
    
    parfor n=1:100
    x01=x1((n-1)*nframes+1:n*nframes);
    s = linspace(0,nframes,nframes*10);
    x2= (sinc_interp(x01.',s)).';
    y01=y1((n-1)*nframes*10+1:n*nframes*10);
    y2 = bandpass(y01,[0.01 0.1],(1/TR)*10); %variable noise
    z01=z1((n-1)*nframes*10+1:n*nframes*10);
    z2 = bandpass(z01,[0.01 0.1],(1/TR)*10); %variable noise
    
    SDx= std(x2);
    SDy= std(y2);
    SDz= std(z2);
    y2 = (y2*SDx./SDy)./5;
    z2 = (z2*SDx./SDz)./5;
    S1=x2+y2*k; % original signal: signal x1+ noisefactor* noise y2
        for i=1:60
            for j=1:60
                    S2= (circshift(x2,i-1))*(1-j/75) + z2*k;
                    R(i,j,n,k)=corr(S1,S2);
            end
        end
    end
    
    mR(:,:,k)=mean(R(:,:,:,k),3);

    %% correlation plot
    figure
    imagesc(mR(:,:,k));
    V=({[0.9 0.9], [0.8,0.8],[0.6 0.6],[0.4 0.4],[0.2 0.2]});
    for l=1:5
    hold on;
    [C,h]=imcontour(mR(:,:,k),V{1 ,l},'w');
    h.LineWidth = 1.5;
    clabel(C,h,'Color','w','FontSize',12);
    end
    colorbar;
    caxis([0.3 1]);
    title(sprintf('Simulated Correlation Noise x%01d',k));
    xlabel('Magnitude Decrease (%)');
    ylabel('Time Delay (msec)')
    yticklabels({'1000','2000','3000','4000','5000','6000'});
    xticks([15 30 45 60])
    xticklabels({'20','40','60','80'});
end
end
