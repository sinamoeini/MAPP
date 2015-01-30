SHELL 	    = /bin/bash
SRC         = src/
PROGRAMS    = main
CC          = mpic++ 
OBJ         = obj/
MAKEFILE    = Makefile
CFLAGS      = -O3 
LIBS        = 
INCLUDES    = 
        
CPP_FILES   = $(wildcard $(SRC)*.cpp)
H_FILES     = $(wildcard $(SRC)*.h)
OBJ_FILES   = $(addprefix $(OBJ),$(notdir $(CPP_FILES:.cpp=.o))) 


$(OBJ)%.o: $(SRC)%.cpp $(MAKEFILE)
	$(CC) -c $(CFLAGS) -o $@ $(INCLUDES) $(LIBS) $<

MAPP:	prep $(OBJ_FILES) $(MAKEFILE) 
	$(CC) $(CFLAGS) $(OBJ_FILES) -o $@ $(INCLUDES) $(LIBS)

clean:  
	rm -rf $(OBJ)
prep:
	@mkdir -p $(OBJ); \
	cd src; \
	rm -rf clock_styles.h; \
	for i in `ls clock_*.h`; do \
	echo \#include \"$$i\" >>clock_styles.h; \
	done; \
	rm -rf command_styles.h; \
	for i in `ls command_*.h`; do \
	echo \#include \"$$i\" >>command_styles.h; \
	done; \
	rm -rf ff_styles.h; \
	for i in `ls ff_*.h`; do \
	echo \#include \"$$i\" >>ff_styles.h; \
	done; \
	rm -rf md_styles.h; \
	for i in `ls md_*.h`; do \
	echo \#include \"$$i\" >>md_styles.h; \
	done; \
	rm -rf min_styles.h; \
	for i in `ls min_*.h`; do \
	echo \#include \"$$i\" >>min_styles.h; \
	done; \
	rm -rf read_styles.h; \
	for i in `ls read_*.h`; do \
	echo \#include \"$$i\" >>read_styles.h; \
	done; \
	rm -rf write_styles.h; \
	for i in `ls write_*.h` ;do \
	echo \#include \"$$i\" >>write_styles.h; \
	done; \
	cd ..
	@echo 'preparing the style files'
