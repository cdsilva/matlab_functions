function data2 = template_align(data, template)

[m, n] = size(data);

data2 = data;

for i=1:m
    temp_data = data(i,:);
    max_ip = -1;
    max_data = 0;
    for j=1:n
        ip = sum(temp_data .* template);
        if ip > max_ip
            max_ip = ip;
            max_data = temp_data;
        end
        temp_data = circshift(temp_data,[0 1]);
    end
    data2(i,:) = max_data;
end
