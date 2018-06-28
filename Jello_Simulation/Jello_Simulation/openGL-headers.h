#if defined(WIN32) || defined(linux)
	#include <GL/gl.h>
	#include <GL/glew.h>
	#include <GL/freeglut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <OpenGL/glu.h>
  #include <GLUT/glut.h>
#endif

