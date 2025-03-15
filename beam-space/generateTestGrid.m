function grid = generateTestGrid(nu, tau, phi, N, M, T, B, fov, Nrf, grid_size)
    df = 1/T;
    dtau = 1/(M*df);
    dnu = 1/(N*T);
    dphi = fov/B/Nrf;

    r1 = (rand(grid_size,1)+1)*dnu/2;
    r1s = rand(grid_size,1);
    r1(r1s<0.5) = -r1(r1s<0.5);

    r2 = (rand(grid_size,1)+1)*dtau/2;
    r2s = rand(grid_size,1);
    r2(r2s<0.5) = -r2(r2s<0.5);

    r3 = (rand(grid_size,1)+1)*dphi/2;
    r3s = rand(grid_size,1);
    r3(r3s<0.5) = -r3(r3s<0.5);

    grid = ones(grid_size,3);
    
    grid(:, 1) = grid(:, 1).*nu + r1;
    grid(:, 2) = grid(:, 2).*tau + r2;
    grid(:, 3) = grid(:, 3).*phi + r3;
    
    grid = grid.';
end