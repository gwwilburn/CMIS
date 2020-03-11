#useful directories
ESLDIR = ../../infernal-1.1.3/easel
HMMERDIR = ../../infernal-1.1.3/hmmer/src
INFERNALDIR = ../../infernal-1.1.3/src
H4DIR = ../../hmmer4/nwo

#variables that I will run
INCLUDES = -I${INFERNALDIR} -I${HMMERDIR} -I${H4DIR} -I${ESLDIR}
IFLAGS   = -I${ESLDIR} -I${HMMERDIR} -I${INFERNALDIR} -I${H4DIR} -I../src
LFLAGS   = -L${ESLDIR} -L${HMMERDIR} -I${INFERNALDIR}  -I${H4DIR} -L../src
LIBS     = -lm
SLIBS    = ${INFERNALDIR}/libinfernal.a ${HMMERDIR}/libhmmer.a ${H4DIR}/libhmmer.a ${ESLDIR}/libeasel.a

#compiler information
CC = gcc
CFLAGS = -Wall -g

all: cmis_scoreset cmis_trace h4_pathalign cmscoreIS cmscoreISh4 cmemit_posterior hmmscore_path

cmis_scoreset: cmis_scoreset.c
	${CC} ${CFLAGS} -c -o cmis_scoreset.o ${IFLAGS} cmis_scoreset.c

cmis_trace: cmis_trace.c
	${CC} ${CFLAGS} -c -o cmis_trace.o ${IFLAGS} cmis_trace.c

h4_pathalign: h4_pathalign.c
	${CC} ${CFLAGS} -c -o h4_pathalign.o ${IFLAGS} h4_pathalign.c

cmscoreIS: cmscoreIS.c
	${CC} ${CFLAGS} -o cmscoreIS ${IFLAGS} ${LFLAGS} ${LIBS} cmscoreIS.c cmis_scoreset.o cmis_trace.o ${SLIBS}

cmscoreISh4: cmscoreISh4.c
	${CC} ${CFLAGS} -o cmscoreISh4 ${IFLAGS} ${LFLAGS} ${LIBS} cmscoreISh4.c cmis_scoreset.o cmis_trace.o h4_pathalign.o ${SLIBS}

cmemit_posterior: cmemit_posterior.c
	${CC} ${CFLAGS} -o cmemit_posterior ${IFLAGS} ${LFLAGS} ${LIBS} cmemit_posterior.c ${SLIBS}

hmmscore_path: hmmscore_path.c
	${CC} ${CFLAGS} -o hmmscore_path ${IFLAGS} ${LFLAGS} ${LIBS} hmmscore_path.c ${SLIBS}
clean:
	rm -f cmscoreIS cmscoreISh4 cmis_scoreset.o cmis_trace.o h4_pathalign.o cmemit_posterior hmmscore_path

.PHONY: all clean
