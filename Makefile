#
# == Paths ==
#
BIN_DIR   := bin
BUILD_DIR := lib
SRC_DIR   := src

#
# == Files ==
#
SRCS := $(wildcard $(SRC_DIR)/*.c)
OBJS  = $(SRCS:.c=.o)
BINS  = $(OBJS:.o=)

#
# == CC Flags ==
#
CC      := mpicc
CFLAGS  := -I${TACC_GSL_INC} -O2 -xHost
LDFLAGS := -L${TACC_GSL_LIB} -lgsl -lgslcblas -qopenmp

#
# == Targets ==
#

all: bins
objs: $(OBJS)
bins: $(BINS)


clean:
	$(RM) $(BUILD_DIR)/*.o $(BIN_DIR)/*

%.o: %.c
	$(CC) -o $(BUILD_DIR)/$(notdir $@) -c $< $(CFLAGS)
%: %.o
	$(CC) -o $(BIN_DIR)/$(basename $(notdir $@)) $(BUILD_DIR)/$(notdir $<) $(LDFLAGS)
