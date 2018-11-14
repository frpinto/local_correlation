function a=checkrow(a)

%this function guarantees that the vector a is in row format and not in
%collumn format.

%Francisco Rodrigues Pinto, Oeiras, 2003


if size (a,2)==1
    a=a';
end