% Tensorise Abilene data by n \times m/T \times T 
function T = mat2ten_day(M, DIV)    
    T= zeros(144,2016/DIV,DIV);
    for i=1:DIV
        s=(i-1)*2016/DIV+1;
        e=i*2016/DIV;
        T(:,:,i) =  M(:,s:e);
    end
end
