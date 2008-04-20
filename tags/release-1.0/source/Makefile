include ../makevars

OBJS = alignment.o alphabet.o amino_acids.o io.o seqio.o similarity.o \
       global_alignment/two_sequences_simple.o global_alignment/two_sequences_linear.o \
       scores/alignment_score.o scores/two_sequences_base.o scores/two_sequences_clusterw.o \
       trees/neighbour_joining_tree.o       
       

$(LIB_NAME): $(OBJS)
	$(AR) lib$(LIB_NAME).a $? 

%.C:

install: $(LIB_NAME)

clean:        
	rm $(OBJS) lib$(LIB_NAME).a
