##
CC = g++
CFLAGS = -O3 -Wall -std=c++11 -I./include

LFLAGS = -pthread

##
SRC = $(wildcard src/*.cpp)
SEARCH = $(wildcard src/searchmain/*.cpp)
INDEX= $(wildcard src/indexmain/*.cpp)
TEST= $(wildcard src/testmain/*.cpp)

OBJS = $(SRC:.cpp=.o)
OBJSEARCH = $(SEARCH:.cpp=.o)
OBJINDEX = $(INDEX:.cpp=.o)
OBJTEST = $(TEST:.cpp=.o)

TARGET1 = nass
TARGET2 = nass-index
TARGET3 = test

##
.suffixes:
	.cpp

%.o:%.cpp
	$(CC) -c $< $(CFLAGS) -o $@

all: $(TARGET1) $(TARGET2) $(TARGET3)

$(TARGET1): $(OBJS) $(OBJSEARCH)
	$(CC) -o $@ $(OBJS) $(OBJSEARCH) $(LFLAGS) 

$(TARGET2): $(OBJS) $(OBJINDEX)
	$(CC) -o $@ $(OBJS) $(OBJINDEX) $(LFLAGS) 

$(TARGET3): $(OBJS) $(OBJTEST)
	$(CC) -o $@ $(OBJS) $(OBJTEST) $(LFLAGS)

clean:
	rm -f *~ */*~ $(OBJS) $(OBJSEARCH) $(OBJINDEX) $(OBJTEST) $(TARGET1) $(TARGET2) $(TARGET3)
