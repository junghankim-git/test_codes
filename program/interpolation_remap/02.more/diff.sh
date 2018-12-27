#!/bin/bash


ncdump result/sr4.2/result.nc > org.txt
ncdump result/result.nc > new.txt
command='diff result/sr4.2/result.nc result/result.nc'
echo $command
diff org.txt new.txt
