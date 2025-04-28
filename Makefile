OBJECTS = config.o grif-replay.o midas-format.o grif-format.o \
          reorder.o default_sort.o test_config.o

CFLAGS  = -g -O2 -fPIC 

grif-replay: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) -rdynamic -lz -ldl -lm -lpthread

midas: midas_module.so

midas_module.so: midas_module.c libmidas.a
	$(CC) $(CFLAGS) -rdynamic -shared -o $@ $^ -lrt -lz -lutil -lnsl -lpthread

.c.o:
	$(CC) -c $(CFLAGS) $<

config.o:       config.c config.h grif-replay.h
grif-format.o:  grif-format.c grif-replay.h grif-format.h midas-format.h
midas-format.o: midas-format.c grif-replay.h midas-format.h
default_sort.o: default_sort.c config.h grif-format.h
grif-replay.o:  grif-replay.c config.h grif-format.h midas-format.h
midas_module.o: midas_module.c grif-format.h config.h midas-format.h midas.h
reorder.o:      reorder.c grif-replay.h midas-format.h

#SOURCES = config.c grif-replay.c midas-format.c grif-format.c \
            reorder.c default_sort.c test_config.c

clean:
	rm -f *.o grif-replay midas_module.so

# if there is a file called "clean", above will fail without this ...
.PHONY: clean
