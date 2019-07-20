function compile_point_set_normals
    mex -v -O point_set_normals.cpp -I/usr/include:/usr/local/include -L/usr/lib:/usr/local/lib -lCGAL -lgmp -lboost_thread
end

