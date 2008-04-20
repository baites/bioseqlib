include ../makevars

EXS = alphabet seqio similarity alignment1 alignment2 global_alignment1 \
      global_alignment2 global_alignment3 score1 score2 score3 score4 score5 trees1 \
      trees2 trees3 trees4 trees5 trees6 trees7 trees8

all: $(EXS)

%.C:

$(EXS): %: %.o
	$(CXX) $@.o -o $@ -L$(LIB_DIR) -l$(LIB_NAME)

clean:        
	rm *.o $(EXS) 
