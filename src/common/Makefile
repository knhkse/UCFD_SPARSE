# Common object Makefile
CC ?= gcc
include $(UCFD_PATH)/Makefile.inc

OBJDIR = ./obj
INVERSE = $(OBJDIR)/inverse.o
CONFIG = ./include/config.h
TYPES = ./include/ucfd_types.h

$(INVERSE) : ./inverse.c $(CONFIG)
	$(CC) $(CFLAGS) -fPIC -I./include -c $< -o $@

common : $(INVERSE) $(CONFIG) $(TYPES)

