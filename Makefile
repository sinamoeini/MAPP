SHELL 	    = /bin/bash
SRC         = src/
PROGRAMS    = main
CC          = mpic++ 
OBJ         = obj/
MAKEFILE    = Makefile
CFLAGS      = -std=c++0x -O3 -fstrict-aliasing -Wstrict-aliasing=2 
LIBS        = 
INCLUDES    = 
        
CPP_FILES   = $(wildcard $(SRC)*.cpp)
H_FILES     = $(wildcard $(SRC)*.h)
OBJ_FILES   = $(addprefix $(OBJ),$(notdir $(CPP_FILES:.cpp=.o))) 


$(OBJ)%.o: $(SRC)%.cpp $(MAKEFILE)
	$(CC) -c $(CFLAGS) -o $@ $(INCLUDES)$<

mapp:	prep $(OBJ_FILES) $(MAKEFILE) 
	$(CC) $(CFLAGS) $(OBJ_FILES) -o $@ $(INCLUDES) $(LIBS)

clean:  
	rm -rf $(OBJ)
prep:
	@mkdir -p $(OBJ); \
	cd src; \
	rm -rf dmd_styles.h; \
	for i in `ls dmd_*.h`; do \
	echo \#include \"$$i\" >>dmd_styles.h; \
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
        rm -rf ls_styles.h; \
        for i in `ls ls_*.h` ;do \
        echo \#include \"$$i\" >>ls_styles.h; \
        done; \
	cd ..
