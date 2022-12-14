#######################################
# HP-Laptop.mk
# Default options for HP-Laptop computer
#######################################
CC=gcc
LIBSLOCAL=-L/lib/x86_64-linux-gnu -llapack -lblas -lm 
INCLUDEBLASLOCAL=-I/usr/include/x86_64-linux-gnu
OPTCLOCAL=-fPIC -march=native