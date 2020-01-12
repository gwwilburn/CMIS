#useful directories
ESLDIR = ../../infernal-1.1.3/easel
HMMERDIR = ../../infernal-1.1.3/hmmer/src
INFERNALDIR = ../../infernal-1.1.3/src

#variables that I will run
INCLUDES = -I${INFERNALDIR} -I${HMMERDIR} -I${ESLDIR}
IFLAGS   = -I${ESLDIR} -I${HMMERDIR} -I${INFERNALDIR} -I../src
LFLAGS   = -L${ESLDIR} -L${HMMERDIR} -I${INFERNALDIR}  -L../src
LIBS     = -lm
SLIBS    = ${INFERNALDIR}/libinfernal.a ${HMMERDIR}/libhmmer.a ${ESLDIR}/libeasel.a

#compiler information
CC = gcc
CFLAGS = -Wall -g

all: cmis_scoreset cmis_trace cmscoreIS

cmis_scoreset: cmis_scoreset.c
	${CC} ${CFLAGS} -c -o cmis_scoreset.o ${IFLAGS} cmis_scoreset.c

cmis_trace: cmis_trace.c
	${CC} ${CFLAGS} -c -o cmis_trace.o ${IFLAGS} cmis_trace.c

cmscoreIS: cmscoreIS.c
	${CC} ${CFLAGS} -o cmscoreIS ${IFLAGS} ${LFLAGS} ${LIBS} cmscoreIS.c cmis_scoreset.o cmis_trace.o ${SLIBS}

clean:
	rm -f cmscoreIS cmis_scoreset.o cmis_trace.o

.PHONY: all clean
