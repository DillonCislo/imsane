function compile_surface_geodesic_pairs

mex -v -O surface_geodesic_pairs.cpp -I/usr/include:/usr/local/include -L/usr/lib:/usr/local/lib -lCGAL -lgmp -lboost_thread -lboost_system

end

