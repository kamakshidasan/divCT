gsx = 256;
gsy = 256;
gsz = 256;

lsx = 128;
lsy = 128;
lsz = 128;
FILE = fopen('Data/foot.raw');
A = fread(FILE,'uint8=>uint8');
B = reshape(A,gsx,gsy,gsz);
S = strcat('parallel',{' '},':::',{' '});
%sd = gs/ls;
 %C = cell((gsx/lsx)^3,1); %assuming each subdomain ios of equal size
 count =1;
%B=F;
for z = 1:lsz:(gsz-lsz+1)
    for y =1:lsy:(gsy-lsy+1)
        for x=1:lsx:(gsx-lsx+1)
            if (x==(gsx-lsx+1)) 
                cx = 0;
            else cx = 1;
            end
            if (y==(gsy-lsy+1)) 
                cy = 0;
            else cy =1;
            end
            if (z==(gsz-lsz+1)) cz = 0;
            else cz =1;
            end
            %fopen'w'
            filestr = (strcat('Data/fo_',num2str(count),'_',num2str(lsy+cy),'_',num2str(lsx+cx),'_',num2str(lsz+cz),'.raw'));
            FILE = fopen(filestr,'w');
            S = strcat(S,'"./serialct',{' '},filestr,{' '}, num2str(lsy+cy),{' '},num2str(lsx+cx),{' '},num2str(lsz+cz),{' '},'f',num2str(count),{' '},num2str(count),'"', {' '}); 
            
            fwrite(FILE,B(y:y+ls-1+cy,x:x+ls-1+cx,z:z+ls-1+cz),'uint8');
            %C{count} = B(y:y+ls-1+cy,x:x+ls-1+cx,z:z+ls-1+cz);
            count = count+1
        end
    end
end
scriptfile = fopen('Data/tree.sh','w');
fprintf(scriptfile,'%s',S{1});
