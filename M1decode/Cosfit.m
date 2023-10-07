function [ f ] = Cosfit( c,x )

%UNTITLED3 Summary of this function goes here

%   Detailed explanation goes here

f=c(1).*cosd(x+c(2))+c(3);
f=f';
end