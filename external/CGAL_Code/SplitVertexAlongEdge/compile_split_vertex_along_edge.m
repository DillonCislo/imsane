function compile_split_vertex_along_edge
    mex -v -O split_vertex_along_edge.cpp -I/usr/include:/usr/local/include -L/usr/lib:/usr/local/lib -lCGAL -lgmp -lboost_thread
end
