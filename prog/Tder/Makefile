OBJS=sample.o VM_vandermonde.o VM_matrix.o
FLAG= -g 

all: sample

sample: $(OBJS)
	$(CXX) -o $@ $(OBJS)

clean: Force
	rm -f $(OBJS) Doxyfile

clean.all: clean
	rm -f sample
	rm -rf html

docs: Force
	./Doxyfile.rb > Doxyfile
	doxygen

Force:
