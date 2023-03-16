if (DOPCC)
SET(clhash_SOURCES 
	clhash.c
)
else()
SET(clhash_SOURCES minim.c)
endif()

ADD_LIBRARY(clhash
${clhash_SOURCES}
)
