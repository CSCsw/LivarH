#include <limits.h>
#include <adolc/adolc.h>

#define NULLLOC UINT_MAX

template <typename T>
class DerivativeInfo {
 public:
  DerivativeInfo() {
    clear();
  }
  ~DerivativeInfo() {
  }

  void clear() {
    r = NULLLOC; x = NULLLOC; y = NULLLOC;
    dx = 0.0; dy = 0.0;
    pxx = 0.0; pxy = 0.0; pyy = 0.0;
    pxxx = 0.0; pxxy = 0.0; pxyy = 0.0; pyyy = 0.0;
  }
  
  void clear(unsigned char opcode) {
    this.opcode = opcode;
    this.clear();
  }

  unsigned char opcode;
  T r,x,y;
  double dx, dy;
  double pxx, pxy, pyy;
  double pxxx, pxxy, pxyy, pyyy;
};
