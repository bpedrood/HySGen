#ifndef LOC_BD_H
#define LOC_BD_H

/////////////////////////////////////////////////
// Bahman Added - to make the code work in Ubuntu 18.04 (Struct is copied from online resources).
#if defined(GLib_GLIBC) || defined(GLib_BSD)
struct __exception {
  int    type;      /* Exception type */
  char*  name;      /* Name of function causing exception */
  double arg1;      /* 1st argument to function */
  double arg2;      /* 2nd argument to function */
  double retval;    /* Function return value */
};
#endif

#endif //LOC_BD_H