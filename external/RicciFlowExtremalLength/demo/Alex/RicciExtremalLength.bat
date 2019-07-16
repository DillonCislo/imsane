set mesh=%1
..\..\bin\RicciFlow.exe -tangent_ricci_extremal_length %mesh%.m %mesh%.uv.m
..\..\bin\SimpleViewer.exe %mesh%.uv.m ..\..\textures\checker_1k.bmp


