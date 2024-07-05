function sides = getHypercubeSides(vertices)
    % vertices: matriks 2^n x n yang mendefinisikan semua titik sudut hyper-kotak
    % sides: array 3D yang menyimpan koordinat titik awal dan akhir dari setiap sisi
    
    num_vertices = size(vertices, 1);  % Jumlah titik sudut
    dim = size(vertices, 2);           % Dimensi ruang
    
    % Inisialisasi list untuk menyimpan sisi-sisi
    sides_list = [];
    
    % Loop untuk setiap pasangan titik sudut
    for i = 1:num_vertices
        for j = i+1:num_vertices
            % Cek jika pasangan titik sudut hanya berbeda di satu dimensi
            if sum(vertices(i,:) ~= vertices(j,:)) == 1
                % Tambahkan sisi ke list
                sides_list = [sides_list; vertices(i,:), vertices(j,:)];
            end
        end
    end
    
    % Konversi list menjadi array 3D
    num_sides = size(sides_list, 1);
    sides = reshape(sides_list', [dim, 2, num_sides]);
    sides = permute(sides, [2, 1, 3]);
end