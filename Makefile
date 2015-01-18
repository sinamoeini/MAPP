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

exec: 	$(OBJ_FILES) $(MAKEFILE)
	$(CC) $(CFLAGS) $(OBJ_FILES) -o $@ $(INCLUDES) $(LIBS)

clean:  
	rm -rf $(OBJ)*.o	
