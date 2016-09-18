CC = mpicc
CFLAGS = -O2 -fopenmp

PROJ = matmultiply_main
OBJS = matmultiply_main.o


matmultiply_main : $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

obj-clean:
	$(RM) *.o

exec-clean:
	$(RM) $(PROJ)

autosave-clean:
	$(RM) *~

clean:
	$(MAKE) obj-clean
	$(MAKE) exec-clean
	$(MAKE) autosave-clean
